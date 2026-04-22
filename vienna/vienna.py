"""
A simple wrapper for RNAfold. Allows folding of RNA sequences to get
secondary structure and energy.
"""

import os
import subprocess
import shutil
import tempfile
import warnings
from dataclasses import dataclass
from typing import List, Tuple

try:
    import RNA

    _HAS_RNA = True
except ImportError as _rna_import_err:
    RNA = None
    _HAS_RNA = False
    warnings.warn(
        f"ViennaRNA Python bindings unavailable ({_rna_import_err}); "
        "falling back to command-line RNAfold/RNAcofold/RNAinverse.",
        stacklevel=2,
    )

# classes #####################################################################


class ViennaException(Exception):
    """
    Exception for Vienna module.
    """


@dataclass(frozen=True, order=True)
class FoldResults:
    """
    Results from calling RNAfold.
    """

    dot_bracket: str
    mfe: float
    ens_defect: float
    bp_probs: List[List[float]]


class InverseResults:
    """
    Results from calling RNAinverse.
    """

    @dataclass(frozen=True, order=True)
    class SeqScore:
        """
        Results from one sequence in the inverse folding.
        """

        seq: str
        score: float

    def __init__(self, seqs: List[str], scores: List[float]) -> None:
        self.seq_scores: List[InverseResults.SeqScore] = [
            self.SeqScore(seq, score) for seq, score in zip(seqs, scores)
        ]

    def __len__(self) -> int:
        return len(self.seq_scores)

    def __iter__(self):
        return iter(self.seq_scores)


# private functions ############################################################
def _get_model_details():
    """
    Get model details matching the command-line options:
    --noLP: no lonely pairs
    -d2: dangles = 2

    Returns:
        RNA.md: Model details object
    """
    md = RNA.md()
    md.noLP = 1  # --noLP: no lonely pairs
    md.dangles = 2  # -d2: dangles = 2
    return md


def _require_cli(tool: str) -> None:
    if shutil.which(tool) is None:
        raise RuntimeError(
            f"{tool} command-line tool not available in PATH and "
            "ViennaRNA Python bindings failed to load."
        )


def _parse_ensemble_diversity(lines: List[str]) -> float:
    for line in lines:
        if "ensemble diversity" in line:
            try:
                return float(line.split()[-1])
            except (ValueError, IndexError):
                return 0.0
    return 0.0


def _fold_cmdline(seq: str, bp_probs: bool) -> "FoldResults":
    _require_cli("RNAfold")
    args = ["RNAfold", "-p", "--noLP", "--noPS", "-d2"]
    if not bp_probs:
        args.append("--noDP")
    with tempfile.TemporaryDirectory() as tmpdir:
        result = subprocess.run(
            args,
            input=seq + "\n",
            capture_output=True,
            text=True,
            cwd=tmpdir,
            check=True,
        )
        lines = result.stdout.split("\n")
        spl1 = lines[1].split()
        structure = spl1[0]
        energy = float(spl1[-1].strip("()"))
        ensemble_diversity = _parse_ensemble_diversity(lines)

        bp_probs_list: List[List[float]] = []
        if bp_probs:
            dot_ps = os.path.join(tmpdir, "dot.ps")
            if os.path.exists(dot_ps):
                with open(dot_ps) as f:
                    for line in f:
                        parts = line.strip().split()
                        # RNAfold writes "i j sqrt(P_ij) ubox" for bp probs
                        if len(parts) == 4 and parts[3] == "ubox":
                            try:
                                i = int(parts[0])
                                j = int(parts[1])
                                p = float(parts[2]) ** 2
                            except ValueError:
                                continue
                            if p > 0.0:
                                bp_probs_list.append([i, j, p])

    return FoldResults(structure, energy, ensemble_diversity, bp_probs_list)


def _cofold_cmdline(seq: str) -> "FoldResults":
    _require_cli("RNAcofold")
    with tempfile.TemporaryDirectory() as tmpdir:
        result = subprocess.run(
            ["RNAcofold", "-p", "--noLP", "--noPS", "-d2"],
            input=seq + "\n",
            capture_output=True,
            text=True,
            cwd=tmpdir,
            check=True,
        )
    lines = result.stdout.split("\n")
    spl1 = lines[1].split()
    structure = spl1[0]
    energy = float(spl1[-1].strip("()"))
    if "&" in seq and "&" not in structure:
        seq1_len = seq.index("&")
        structure = structure[:seq1_len] + "&" + structure[seq1_len:]
    return FoldResults(structure, energy, _parse_ensemble_diversity(lines), [])


def _score_cmdline(seq: str) -> float:
    _require_cli("RNAfold")
    result = subprocess.run(
        ["RNAfold", "--noLP", "--noPS", "-d2"],
        input=seq + "\n",
        capture_output=True,
        text=True,
        check=True,
    )
    lines = result.stdout.split("\n")
    return float(lines[1].split()[-1].strip("()"))


# public functions #############################################################
def fold(seq: str, bp_probs: bool = False) -> FoldResults:
    """
    Fold an RNA sequence using RNAfold.

    Args:
        seq (str): The RNA sequence to fold.
        bp_probs (bool): Generate base pair probabilities? (default: False)

    Returns:
        FoldResults: Results from RNAfold
    """
    if len(seq) == 0:
        raise ValueError("Must supply a sequence longer than 0")

    if not _HAS_RNA:
        return _fold_cmdline(seq, bp_probs)

    md = _get_model_details()

    # Create fold compound with model details
    fc = RNA.fold_compound(seq, md)

    # Compute MFE structure and energy
    structure, energy = fc.mfe()

    # Calculate partition function for ensemble diversity and base pair probabilities
    fc.pf()

    # Calculate ensemble diversity (mean pairwise Hamming distance between structures)
    # This is what RNAfold outputs as "ensemble diversity"
    # Note: This is different from ensemble_defect, which is the distance from MFE structure
    try:
        ensemble_diversity = fc.mean_bp_distance()
    except (AttributeError, TypeError):
        # Fallback: use 0.0 if method not available
        ensemble_diversity = 0.0

    bp_probs_list = []
    if bp_probs:
        # Get base pair probability matrix
        bpp_matrix = fc.bpp()
        n = len(seq)
        # bpp_matrix uses 0-based indexing: bpp[i][j] corresponds to positions i+1, j+1
        for i in range(n):
            for j in range(i + 1, n):
                prob = bpp_matrix[i][j]
                if prob > 0.0:
                    bp_probs_list.append(
                        [i + 1, j + 1, prob]
                    )  # +1 for 1-indexed output

    return FoldResults(structure, energy, ensemble_diversity, bp_probs_list)


def cofold(seq: str) -> FoldResults:
    """
    Cofold two RNA sequences to get their combined secondary structure and energy.

    Args:
        seq (str): Sequences to fold separated by a '&'

    Returns:
        FoldResults: Results from RNAcofold
    """
    if len(seq) == 0:
        raise ValueError("Must supply a sequence longer than 0")

    if "&" not in seq:
        raise ValueError("Cofold requires sequences separated by '&'")

    if not _HAS_RNA:
        return _cofold_cmdline(seq)

    seq1, seq2 = seq.split("&", 1)

    md = _get_model_details()

    # Create fold compound with the combined sequence (RNA.fold_compound handles "&" separator)
    fc = RNA.fold_compound(seq, md)

    # Compute MFE structure and energy
    structure, energy = fc.mfe()

    # Insert "&" separator in structure at position corresponding to seq1 length
    # This matches command-line RNAcofold output format
    structure = structure[: len(seq1)] + "&" + structure[len(seq1) :]

    # Calculate partition function for ensemble diversity
    # Note: This can be slow for long sequences, but we include it for consistency
    fc.pf()

    # Calculate ensemble diversity (mean pairwise Hamming distance between structures)
    try:
        ensemble_diversity = fc.mean_bp_distance()
    except (AttributeError, TypeError):
        ensemble_diversity = 0.0

    return FoldResults(structure, energy, ensemble_diversity, [])


def inverse_fold(secstruct: str, constraint: str, n_sol: int = 100) -> InverseResults:
    """
    Generates sequences that match a secondary structure with sequence constraint.

    Args:
        secstruct (str): Secondary structure in dot bracket notation
        constraint (str): Sequence constraints
        n_sol (int): Number of solutions to return (default: 100)

    Returns:
        InverseResults: Results from RNAinverse
    """
    # Check if RNAinverse is available
    if shutil.which("RNAinverse") is None:
        raise RuntimeError("RNAinverse command-line tool not available in PATH")

    md = _get_model_details() if _HAS_RNA else None
    seqs = []
    scores = []

    # Call RNAinverse with -R option to get multiple solutions
    # Format: structure on first line, constraint on second line
    input_data = f"{secstruct}\n{constraint}\n"
    try:
        result = subprocess.run(
            ["RNAinverse", f"-R{n_sol}"],
            input=input_data,
            capture_output=True,
            text=True,
            check=True,
        )

        # Parse output: each line is "SEQUENCE    DISTANCE"
        # Include all sequences returned (distance indicates how close they are to target)
        lines = result.stdout.strip().split("\n")
        for line in lines:
            if not line.strip():
                continue

            # Split by whitespace - first part is sequence, second is distance (if present)
            parts = line.strip().split()
            if parts:
                seq = parts[0]
                if seq and seq not in seqs:  # Avoid duplicates
                    seqs.append(seq)

                    # Compute score by folding the sequence
                    if _HAS_RNA:
                        fc = RNA.fold_compound(seq, md)
                        _, mfe = fc.mfe()
                    else:
                        mfe = _score_cmdline(seq)
                    scores.append(mfe)

                    # Stop if we have enough solutions
                    if len(seqs) >= n_sol:
                        break

    except subprocess.CalledProcessError as e:
        # If RNAinverse fails, raise an error
        raise RuntimeError(f"RNAinverse failed: {e.stderr}") from e
    except Exception as e:
        raise RuntimeError(f"Error running RNAinverse: {e}") from e

    return InverseResults(seqs, scores)


def folded_structure(seq: str) -> str:
    """
    Get just the folded structure

    Args:
        seq (str): RNA sequence
    """
    fold_result = fold(seq)
    return fold_result.dot_bracket


def does_sequence_fold_to(seq: str, target_structure: str) -> bool:
    """
    Determines if a sequence folds into a specific secondary structure.

    Args:
        seq (str): RNA sequence
        target_structure (str): Target secondary structure in dot bracket notation

    Returns:
        bool: True if the sequence folds into the target structure, False otherwise
    """
    actual_structure = folded_structure(seq)
    return actual_structure == target_structure
