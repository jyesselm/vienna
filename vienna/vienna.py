import os
import subprocess
import shutil
from dataclasses import dataclass
from typing import List, Optional

RNA_FOLD_EXISTS = False


class ViennaException(Exception):
    pass


@dataclass(frozen=True, order=True)
class FoldResults(object):
    dot_bracket: str
    mfe: float
    ens_defect: float
    bp_probs: List[List[float]]


class InverseResults(object):
    class SeqScore(object):
        def __init__(self, seq, score):
            self.seq = seq
            self.score = score

    def __init__(self, seqs, scores):
        self.seq_scores = []
        for seq, score in zip(seqs, scores):
            self.seq_scores.append(self.SeqScore(seq, score))

    def __len__(self):
        return len(self.seq_scores)

    def __iter__(self):
        return self.seq_scores.__iter__()


def fold(seq: str, bp_probs=False) -> FoldResults:
    global RNA_FOLD_EXISTS
    if not RNA_FOLD_EXISTS:
        if shutil.which("RNAfold") is None:
            raise ViennaException("RNAfold is not in the path!")
        else:
            RNA_FOLD_EXISTS = True

    if len(seq) == 0:
        raise ValueError("must supply a sequence longer then 0")
    if bp_probs:
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
    spl1 = lines[1].split()
    spl2 = lines[-2].split()
    ensemble_diversity = float(spl2[-1])
    structure = spl1[0]
    energy = float(lines[1].split("(")[-1][:-1])
    bp_probs_list = []
    if bp_probs:
        f = open("dot.ps")
        lines = f.readlines()
        f.close()
        for l in lines:
            spl = l.split()
            if len(spl) != 4:
                continue
            if spl[3] != "ubox":
                continue
            bp_probs_list.append([int(spl[0]), int(spl[1]), float(spl[2])])

    results = FoldResults(structure, energy, ensemble_diversity, bp_probs_list)
    try:
        os.remove("rna.ps")
        os.remove("dot.ps")
    except:
        pass
    return results


def folded_structure(seq: str) -> str:
    r = fold(seq)
    return r.dot_bracket


def cofold(seq):
    global RNA_FOLD_EXISTS
    if not RNA_FOLD_EXISTS:
        if shutil.which("RNAfold") is None:
            raise ViennaException("RNAfold is not in the path!")
        else:
            RNA_FOLD_EXISTS = True

    if len(seq) == 0:
        raise ValueError("must supply a sequence longer then 0")
    os.system('echo "' + seq + '" | ' + "RNAcofold -p > rnafold_dump")
    f = open("rnafold_dump")
    lines = f.readlines()
    f.close()
    try:
        last_line = lines.pop()
    except:
        return None
    spl = last_line.split()
    ensemble_prob = float(spl[6][:-1])
    spl = lines[1].split()
    spl2 = lines[1].split("(")
    structure = spl[0]
    energy = float(spl2[-1][:-2].rstrip())
    results = FoldResults(structure, energy, ensemble_prob)
    try:
        os.remove("rnafold_dump")
        os.remove("rna.ps")
        os.remove("dot.ps")
    except:
        pass
    return results


def inverse(ss, constraint, n_sol=100, discard_misfolds=True):
    global RNA_FOLD_EXISTS
    if not RNA_FOLD_EXISTS:
        if shutil.which("RNAfold") is None:
            raise ViennaException("RNAfold is not in the path!")
        else:
            RNA_FOLD_EXISTS = True

    f = open("seqsecstruct.txt", "w")
    f.writelines([ss + "\n", constraint])
    f.close()

    try:
        output = subprocess.check_output(
            "RNAinverse -Fmp -f 0.5 -d2 -R{} < seqsecstruct.txt".format(n_sol),
            shell=True,
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            "command '{}' return with error (code {}): {}".format(
                e.cmd, e.returncode, e.output
            )
        )
    lines = output.decode("utf-8").split("\n")
    seqs = []
    scores = []
    for e in lines:
        print(e)
        spl = e.split()
        if len(spl) != 2:
            continue
        seqs.append(spl[0])
        scores.append(int(spl[1]))
    return InverseResults(seqs, scores)
