"""
A simple wrapper for RNAfold. Allows folding of RNA sequences to get
secondary structure and energy.
"""

import re
import os
import subprocess
import shutil
from dataclasses import dataclass
from typing import List, Tuple

# classes #####################################################################


@dataclass(order=True)
class Globals:
    """
    Global variables for module
    """

    rna_fold_exists: bool = False
    rna_cofold_exists: bool = False
    rna_inverse_exists: bool = False
    version: str = ""


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


# module vars ##################################################################

globs = Globals()


# private functions ############################################################
def _get_fold_results(lines: List[str]) -> Tuple[float, str, float]:
    """
    Get results from RNAfold output.

    Args:
        lines (List[str]): lines from RNAfold output

    Returns:
        Tuple[float, str, float]: ensemble defect, structure, and energy
    """
    spl1 = lines[1].split()
    spl2 = lines[-2].split()
    try:
        ensemble_diversity = float(spl2[-1])
    except ValueError:
        ensemble_diversity = 0.0
    mfe = float(spl1[-1].strip("()"))
    structure = spl1[0]
    return ensemble_diversity, structure, mfe


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
    if not globs.rna_fold_exists:
        if shutil.which("RNAfold") is None:
            raise ViennaException("RNAfold is not in the path!")
        globs.rna_fold_exists = True

    if globs.version == "":
        output = subprocess.check_output("RNAfold --version", shell=True)
        spl = output.decode("utf-8").split()
        globs.version = spl[1]

    if len(seq) == 0:
        raise ValueError("Must supply a sequence longer than 0")

    ver_spl = globs.version.split(".")
    cmd = (
        f'echo "{seq}" | RNAfold -p --noLP --noPS -d2'
        if bp_probs or int(ver_spl[1]) < 5
        else f'echo "{seq}" | RNAfold -p --noLP --noDP --noPS -d2'
    )

    output = subprocess.check_output(cmd, shell=True)
    lines = output.decode("utf-8").split("\n")
    ens_defect, structure, energy = _get_fold_results(lines)
    bp_probs_list = []

    if bp_probs:
        with open("dot.ps", "r", encoding="UTF-8") as fhandler:
            lines = fhandler.readlines()
        for line in lines:
            spl = line.split()
            if len(spl) != 4:
                continue
            if spl[3] != "ubox":
                continue
            bp_probs_list.append([int(spl[0]), int(spl[1]), float(spl[2])])
        try:
            os.remove("dot.ps")
        except FileNotFoundError:
            pass

    return FoldResults(structure, energy, ens_defect, bp_probs_list)


def cofold(seq: str) -> FoldResults:
    """
    Cofold two RNA sequences to get their combined secondary structure and energy.

    Args:
        seq (str): Sequences to fold separated by a '&'

    Returns:
        FoldResults: Results from RNAcofold
    """
    if not globs.rna_cofold_exists:
        if shutil.which("RNAcofold") is None:
            raise ViennaException("RNAcofold is not in the path!")
        globs.rna_cofold_exists = True

    if len(seq) == 0:
        raise ValueError("Must supply a sequence longer than 0")

    output = subprocess.check_output(
        f'echo "{seq}" | RNAcofold -p --noLP --noPS -d2', shell=True
    )
    lines = output.decode("utf-8").split("\n")
    ens_defect, structure, energy = _get_fold_results(lines)
    return FoldResults(structure, energy, ens_defect, [])


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
    if not globs.rna_inverse_exists:
        if shutil.which("RNAinverse") is None:
            raise ViennaException("RNAinverse is not in the path!")
        globs.rna_inverse_exists = True

    with open("inverse.in", "w", encoding="utf-8") as fhandler:
        fhandler.write(f"{secstruct}\n{constraint}\n")

    output = subprocess.check_output(
        f"RNAinverse -Fmp -f 0.5 -d2 -R{n_sol} < inverse.in", shell=True
    )
    lines = output.decode("utf-8").split("\n")
    seqs = []
    scores = []
    for line in lines:
        spl = line.split()
        if len(spl) != 2:
            continue
        seqs.append(spl[0])
        scores.append(float(spl[1]))

    try:
        os.remove("inverse.in")
        os.remove("dot.ps")
    except FileNotFoundError:
        pass  # ignore

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
    return folded_structure(seq, target_structure)
