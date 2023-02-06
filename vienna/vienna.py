"""
A simple wrapper for RNAfold. Allows folding of RNA sequences to get
secondary structure and energy.
"""

import os
import subprocess
import shutil
from dataclasses import dataclass
from typing import List


# classes #####################################################################

# TODO add check for version if below a set version do not use -noDP flag or it will error out
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

    def __init__(self, seqs, scores):
        self.seq_scores = []
        for seq, score in zip(seqs, scores):
            self.seq_scores.append(self.SeqScore(seq, score))

    def __len__(self):
        return len(self.seq_scores)

    def __iter__(self):
        return self.seq_scores.__iter__()


# module vars ##################################################################

globs = Globals()


# private functions ############################################################
def _get_fold_results(lines):
    """
    Get results from RNAfold output.
    :param lines: lines from RNAfold output
    :return: ens_defect, structure, energy
    """
    spl1 = lines[1].split()
    spl2 = lines[-2].split()
    ensemble_diversity = float(spl2[-1])
    structure = spl1[0]
    energy = float(lines[1].split("(")[-1][:-1])
    return ensemble_diversity, structure, energy


# functions ####################################################################


def fold(seq: str, bp_probs=False) -> FoldResults:
    """
    Fold a sequence using RNAfold.
    :param seq: sequence to fold.
    :param bp_probs: generate base pair probabilities?
    :return: Results from RNAfold [FoldResults] object
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
        raise ValueError("must supply a sequence longer then 0")
    ver_spl = globs.version.split(".")
    if bp_probs or int(ver_spl[1]) < 5:
        output = subprocess.check_output(
            'echo "' + str(seq) + '" | ' + "RNAfold -p --noLP --noPS -d2",
            shell=True,
        )
    else:
        output = subprocess.check_output(
            'echo "' + str(seq) + '" | ' + "RNAfold -p --noLP --noDP --noPS -d2",
            shell=True,
        )
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


def folded_structure(seq: str) -> str:
    """
    Fold a sequence using RNAfold and return the structure.
    :param seq: sequence to fold
    :return: structure in dot-bracket notation
    """
    return fold(seq).dot_bracket


def does_sequence_fold_to(seq: str, struct: str) -> bool:
    """
    Check if a sequence folds to a given structure.
    :param seq: sequence to check
    :param struct: structure to check
    :return: True if sequence folds to structure, False otherwise
    """
    return folded_structure(seq) == struct


def cofold(seq: str) -> FoldResults:
    """
    fold sequences using RNAcofold.
    :param seq: sequences to fold seperated by a '&'
    :return: Results from RNAfold [FoldResults] object
    """
    if not globs.rna_cofold_exists:
        if shutil.which("RNAcofold") is None:
            raise ViennaException("RNAcofold is not in the path!")
        globs.rna_cofold_exists = True
    if len(seq) == 0:
        raise ValueError("must supply a sequence longer then 0")
    output = subprocess.check_output(
        'echo "' + str(seq) + '" | ' + "RNAcofold -p --noLP --noPS -d2",
        shell=True,
    )
    lines = output.decode("utf-8").split("\n")
    ens_defect, structure, energy = _get_fold_results(lines)
    return FoldResults(structure, energy, ens_defect, [])


def inverse_fold(secstruct, constraint, n_sol=100) -> InverseResults:
    """
    Generates sequences that match a secondary structure with sequence constraint
    For example:
    secstruct:  (((.(((....))).)))
    constraint: NNNgNNNNNNNNNNaNNN
    :param secstruct: secondary structure in dot bracket notation
    :param constraint: sequence constraints
    :param n_sol: number of solutions to return
    :return: InverseResults object
    """
    if not globs.rna_inverse_exists:
        if shutil.which("RNAinverse") is None:
            raise ViennaException("RNAinverse is not in the path!")
        globs.rna_inverse_exists = True
    with open("inverse.in", "w", encoding="utf-8") as fhandler:
        fhandler.write(secstruct + "\n" + constraint + "\n")
    output = subprocess.check_output(
        f"RNAinverse -Fmp -f 0.5 -d2 -R{n_sol} < inverse.in",
        shell=True,
    )
    lines = output.decode("utf-8").split("\n")
    seqs = []
    scores = []
    for line in lines:
        spl = line.split()
        if len(spl) != 2:
            continue
        seqs.append(spl[0])
        scores.append(int(spl[1]))
    try:
        os.remove("inverse.in")
        os.remove("dot.ps")
    except FileNotFoundError:
        pass  # ignore
    return InverseResults(seqs, scores)
