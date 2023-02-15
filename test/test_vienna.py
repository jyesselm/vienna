"""
Tests for `vienna` vienna.py module.
"""
from vienna import fold, folded_structure, cofold, inverse_fold


def test_fold():
    """
    Test the fold function
    """
    r = fold("GGGGAAAACCCC")
    assert r.dot_bracket == "((((....))))"


def test_fold_longer():
    """
    Test the fold function
    """
    r = fold("GGGGAAAACCCC" * 20)
    assert r.mfe < -190

def test_fold_bp_prob():
    """
    Test the bp_probs in the fold function
    """
    r = fold("GGGGAAAACCCC", bp_probs=True)
    assert len(r.bp_probs) == 14


def test_folded_structure():
    """
    Test the folded structure function
    """
    struct = folded_structure("GGGGAAAACCCC")
    assert struct == "((((....))))"


def test_cofold():
    """
    Test the cofold function
    """
    r = cofold("GGGG&AAACCCC")
    assert r.dot_bracket == "((((&...))))"


def test_cofold_longer():
    seq = "G" * 20 + "&" + "C" * 20
    r = cofold(seq)


def test_inverse_fold():
    """
    Test the inverse fold function
    """
    r = inverse_fold("(((.(((....))).)))", "NNNgNNNNNNNNNNaNNN", n_sol=5)
    assert len(r) == 5
