"""
Tests for `vienna` module.
"""
import pytest
from vienna import vienna

def test():
    r1 = vienna.fold('GGGGAAAACCCC')
    assert r1.dot_bracket == '((((....))))'
