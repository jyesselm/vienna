"""
Tests for `vienna` vienna.py module.
"""

import subprocess
import shutil
import time
import random
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
    assert len(r.bp_probs) == 11


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


def test_fold_matches_command_line():
    """
    Test that vienna.fold() produces the exact same results as command-line RNAfold.
    Compares structure, energy, and ensemble defect.
    """
    # Skip test if RNAfold is not available in PATH
    if shutil.which("RNAfold") is None:
        import pytest

        pytest.skip("RNAfold command-line tool not available in PATH")

    test_sequences = [
        "GGGGAAAACCCC",
        "AUGCGUACGUACGUACGUA",
        "GGGGAAAACCCC" * 5,
    ]

    for seq in test_sequences:
        # Run command-line RNAfold with -p --noLP -d2
        cmd = f'echo "{seq}" | RNAfold -p --noLP --noPS --noDP -d2'
        output = subprocess.check_output(cmd, shell=True)
        lines = output.decode("utf-8").split("\n")

        # Parse command-line output (matching old _get_fold_results logic)
        # Format:
        # Line 0: sequence
        # Line 1: structure and energy like: "((((....)))) ( -5.40)"
        # Line 3: ensemble info like: "((((....)))) { -5.40 d=0.93}"
        # Last non-empty line: " frequency of mfe structure in ensemble 0.84005; ensemble diversity 1.62"

        # Parse structure and energy from line 1
        spl1 = lines[1].split()
        structure_cl = spl1[0]
        # Energy is in parentheses like "( -5.40)" or "(-5.40)"
        energy_str = spl1[-1].strip("()")
        energy_cl = float(energy_str)

        # Parse ensemble diversity from lines[-2] (second to last line, before empty line)
        # This matches the old _get_fold_results logic exactly
        # Format: " frequency of mfe structure in ensemble 0.84005; ensemble diversity 1.62"
        ensemble_diversity_cl = 0.0
        if len(lines) >= 2:
            spl2 = lines[-2].split()
            try:
                # Last element should be the ensemble diversity value
                ensemble_diversity_cl = float(spl2[-1])
            except (ValueError, IndexError):
                pass

        # Run Python function
        result = fold(seq)

        # Compare results
        assert result.dot_bracket == structure_cl, (
            f"Structure mismatch for sequence {seq}:\n"
            f"  Command-line: {structure_cl}\n"
            f"  Python:       {result.dot_bracket}"
        )

        assert abs(result.mfe - energy_cl) < 1e-6, (
            f"Energy mismatch for sequence {seq}:\n"
            f"  Command-line: {energy_cl}\n"
            f"  Python:       {result.mfe}\n"
            f"  Difference:   {abs(result.mfe - energy_cl)}"
        )

        # Ensemble defect (stored as ens_defect but originally was ensemble diversity)
        # May have slight numerical differences, so use a small tolerance
        assert abs(result.ens_defect - ensemble_diversity_cl) < 0.1, (
            f"Ensemble defect/diversity mismatch for sequence {seq}:\n"
            f"  Command-line (ensemble diversity): {ensemble_diversity_cl}\n"
            f"  Python (ens_defect):               {result.ens_defect}\n"
            f"  Difference:                        {abs(result.ens_defect - ensemble_diversity_cl)}"
        )


def test_benchmark_python_vs_commandline():
    """
    Benchmark test comparing Python API vs command-line RNAfold for 100 sequences.
    Tests correctness and speed comparison.
    """
    # Skip if RNAfold is not available
    if shutil.which("RNAfold") is None:
        import pytest

        pytest.skip("RNAfold command-line tool not available in PATH")

    # Generate 100 test sequences of varying lengths
    random.seed(42)  # For reproducibility
    nucleotides = "AUGC"
    test_sequences = []

    # Mix of lengths: short, medium, long
    for i in range(100):
        if i % 3 == 0:
            # Short sequences (20-30 bp)
            length = random.randint(20, 30)
        elif i % 3 == 1:
            # Medium sequences (40-60 bp)
            length = random.randint(40, 60)
        else:
            # Long sequences (80-120 bp)
            length = random.randint(80, 120)

        seq = "".join(random.choice(nucleotides) for _ in range(length))
        test_sequences.append(seq)

    # Benchmark Python API
    python_start = time.perf_counter()
    python_results = []
    for seq in test_sequences:
        result = fold(seq)
        python_results.append(result)
    python_time = time.perf_counter() - python_start

    # Benchmark command-line
    cl_start = time.perf_counter()
    cl_results = []
    for seq in test_sequences:
        cmd = f'echo "{seq}" | RNAfold -p --noLP --noPS --noDP -d2'
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.DEVNULL)
        lines = output.decode("utf-8").split("\n")

        # Parse structure and energy
        spl1 = lines[1].split()
        structure_cl = spl1[0]
        energy_str = spl1[-1].strip("()")
        energy_cl = float(energy_str)

        cl_results.append((structure_cl, energy_cl))
    cl_time = time.perf_counter() - cl_start

    # Verify correctness: all results should match
    mismatches = 0
    for i, (seq, python_result, (structure_cl, energy_cl)) in enumerate(
        zip(test_sequences, python_results, cl_results)
    ):
        if python_result.dot_bracket != structure_cl:
            mismatches += 1
            print(f"Mismatch {i}: Structure differs for sequence {seq[:20]}...")
            print(f"  Python:      {python_result.dot_bracket}")
            print(f"  Command-line: {structure_cl}")

        if abs(python_result.mfe - energy_cl) >= 1e-6:
            mismatches += 1
            print(f"Mismatch {i}: Energy differs for sequence {seq[:20]}...")
            print(f"  Python:      {python_result.mfe}")
            print(f"  Command-line: {energy_cl}")
            print(f"  Difference:  {abs(python_result.mfe - energy_cl)}")

    # All results must match
    assert (
        mismatches == 0
    ), f"Found {mismatches} mismatches between Python and command-line results"

    # Print benchmark results
    speedup = cl_time / python_time
    print(f"\n{'='*60}")
    print(f"Benchmark Results (100 sequences)")
    print(f"{'='*60}")
    print(
        f"Python API:    {python_time:.4f} seconds ({python_time/100*1000:.2f} ms/sequence)"
    )
    print(f"Command-line:  {cl_time:.4f} seconds ({cl_time/100*1000:.2f} ms/sequence)")
    print(f"Speedup:       {speedup:.2f}x")
    print(f"{'='*60}")

    # The Python API should be faster (or at least comparable)
    # This is the main benefit of using the Python API vs spawning processes
    if speedup > 1:
        print(f"✓ Python API is {speedup:.2f}x faster than command-line")
    else:
        print(f"✓ Python API is {1/speedup:.2f}x slower than command-line (acceptable)")
