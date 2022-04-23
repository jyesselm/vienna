"""
Tests for `vienna` module.
"""
import pytest
from vienna import vienna

def test():
    r1 = vienna.fold('GGGGAAAACCCC')
    assert r1.dot_bracket == '((((....))))'

def test_bp_prob():
    r1 = vienna.fold('GGGGAAAACCCC', bp_probs=True)
    assert len(r1.bp_probs) == 14

