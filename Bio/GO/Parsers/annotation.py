#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
Parser for Gene Ontology annotation files

Usage example::

    import Bio.GO.Parsers.oboparser as obo
    import Bio.GO.Parsers.annotation as annotation

    parser = obo.Parser(open("gene_ontology.1_2.obo"))
    ontology = parser.parse()

    parser = annotation.Parser(open("gene_associaiton.sgd"))
    annotations = parser.parse(ontology)
"""

__author__  = "Tamas Nepusz"
__email__   = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__version__ = "0.1"


__all__ = ["ParseError", "Annotation", "Parser"]

from Bio.GO.ontology import Ontology


class ParseError(Exception):
    """Exception thrown when a parsing error occurred"""

    def __init__(self, msg, lineno = 1):
        Exception.__init__("%s near line %d" % (msg, lineno))
        self.lineno = lineno


class Annotation(object):
    """Class representing an entry from an annotation file.

    We support both GAF 1.0 and 2.0 annotation file formats. These files
    are simple flat files with 15 (GAF 1.0) or 17 (GAF 2.0) tab-delimited
    fields, some of which are required for all entries. Lines starting
    with an exclamation mark are silently ignored.

    An annotation consists of the following fields:

    """

    __slots__ = ["name", "tags"]

    def __init__(self, name, tags=None):
        """Creates a new stanza with the given name and the given
        tags (which must be a dict)"""
        self.name = name
        if tags:
            self.tags = dict(tags)
        else:
            self.tags = dict()

    def __repr__(self):
        """Returns a Python representation of this object"""
        return "%s(%r, %r)" % (self.__class__.__name__, \
                self.name, self.tags)

    def add_tag(self, name, value):
        """Adds a tag-value pair to this stanza. If the tag name already
        exists, the value will be appended to the value list of that tag.
        """
        try:
            self.tags[name].append(value)
        except KeyError:
            self.tags[name] = [value]

    def to_term(self, term_factory=GOTerm):
        """Converts this stanza to an instance of the given term class.

        The stanza must contain a term; in other words, the stanza name must
        be equal to ``"Term"``. This is not checked in this method.

        :Parameters:
        - `term_factory`: the term factory to be used. It is safe to leave
          it at its default value unless you want to use a custom term class.
          The signature of this factory function should be identical to that
          of the constructor of `GOTerm`.
        """
        identifier = str(self.tags["id"][0])
        name = str(self.tags.get("name", [identifier])[0])
        aliases = [str(go_id) for go_id in self.tags.get("alt_id", [])]
        return term_factory(identifier, name, aliases, self.tags)


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
    parser as if it were a list. The iterator yields `Annotation`
    objects. If you are not interested in the individual `Annotation`s
    and you only need an `Ontology` object, call the `parse()`
    method which will construct it for you.
    """

    def __init__(self, fp):
        """Creates an annotation parser that reads the given file-like object.
        """
        self._line_buffer = []

        if isinstance(fp, (str, unicode)):
            fp = open(fp)
        self.fp = fp
        self.lineno = 0
        self._read_headers()

    def _lines(self):
        """Iterates over the lines of the file, removing
        comments and trailing newlines and merging multi-line
        tag-value pairs into a single line"""
        while self._line_buffer:
            yield self._line_buffer.pop(-1)

        while True:
            self.lineno += 1
            line = self.fp.readline()
            if not line:
                break

            line = line.strip()
            if not line:
                # This is an empty line, so it can be ignored
                continue

            yield line

    def _read_headers(self):
        """Reads the headers from the annotation file"""
        self.headers = {}
        for line in self._lines():
            if line[0] != '!':
                # We have reached the end of headers
                self._line_buffer.append(line)
                return
            key, value = [part.strip() for part in line[1:].split(":", 1)]
            try:
                self.headers[key].append(value)
            except KeyError:
                self.headers[key] = [value]

    def annotations(self):
        """Iterates over the annotations in this annotation file,
        yielding an `Annotation` object for each annotation."""
        annotation = None

        for line in self._lines():
            if line[0] == '!':
                continue
            parts = line.strip().split()
            # TODO
            yield Annotation()

    def __iter__(self):
        return self.annotations()

    def parse(self, ontology):
        """Parses the file handle given during construction time and
        returns an appropriately constructed `Annotations` instance.
        
        :Parameters:
            - `ontology`: the ontology being used to map term IDs to
              term names
        """
        for stanza in self:
            # TODO
            pass
