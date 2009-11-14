#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
A parser for the OBO format.

"""

__author__ = 'Chris Lasher'
__email__ = 'chris DOT lasher <AT> gmail DOT com'

from Bio.GO import ontology

from pyparsing import (
        Word,
        alphas,
        nums,
        SkipTo,
        Suppress,
        restOfLine,
        Optional,
        Group,
        Literal,
        lineEnd,
        ZeroOrMore
        )

# comments begin with an exclamation point; anything after them needs to
# be parsed as a comment and otherwise ignored
comment = Suppress('!') + SkipTo(lineEnd)
comment.setParseAction(lambda t: t[0].strip())
# the parsing portion of a line should end either at a newline, or at
# the appearance of a comment
end = comment | lineEnd
# value strings can be broken up over multiple lines if escaped by a
# backslash character immediately before the end of a line (newline)
end_escape = Literal('\\') + lineEnd
continuation = Suppress(end_escape) + ZeroOrMore(SkipTo(end_escape)) +\
        SkipTo(end)
tag = Word(alphas + '_-')
tag.setParseAction(lambda tokens: ' '.join(tokens))
value = (SkipTo(continuation) + continuation) | SkipTo(end)
# we want the value returned as a single string, rather than a disjoint
# list of strings
value.setParseAction(lambda tokens: ''.join((t.lstrip() for t in
    tokens)))
tag_value_pair = tag('tag') + Suppress(':') + value('value') + \
        Optional(comment('comment'))
