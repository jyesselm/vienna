# vienna

[![PYPI package](https://badge.fury.io/py/vienna.png)](http://badge.fury.io/py/vienna)
[![Build status](https://travis-ci.org/jyesselm/vienna.png?branch=main)](https://travis-ci.org/jyesselm/vienna)
[![linting: pylint](https://img.shields.io/badge/linting-pylint-yellowgreen)](https://github.com/PyCQA/pylint)
[![formatting: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

a short python wrapper for vienna
tools (https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/install.html) Note I did not develop this
software just built this short wrapper.

## Installation

### Using Conda (Recommended)

The easiest way to install is using conda with the provided environment file:

```shell
# Clone the repository
git clone https://github.com/jyesselm/vienna
cd vienna

# Create and activate the conda environment (includes ViennaRNA)
conda env create -f environment.yml
conda activate vienna
```

### Using pip

Note: This is just a wrapper, so you must install ViennaRNA first.

**Install ViennaRNA:**
```shell
# Using conda (recommended)
conda install -c bioconda viennarna

# Or install from source: https://www.tbi.univie.ac.at/RNA/
```

**Install the Python package:**
```shell
# From PyPI
pip install vienna

# Or from source
git clone https://github.com/jyesselm/vienna
cd vienna
pip install .
```

## Usage

### Basic Folding

The simplest way to fold an RNA sequence:

```python
import vienna

# Fold a simple RNA sequence
result = vienna.fold('GGGGAAAACCCC')
print(result)
# FoldResults(dot_bracket='((((....))))', mfe=-5.4, ens_defect=1.62)

# Access individual properties
print(f"Structure: {result.dot_bracket}")
print(f"Minimum Free Energy: {result.mfe} kcal/mol")
print(f"Ensemble Diversity: {result.ens_defect}")
```

### Getting Just the Structure

If you only need the secondary structure:

```python
structure = vienna.folded_structure('GGGGAAAACCCC')
print(structure)  # '((((....))))'
```

### Folding with Base Pair Probabilities

Calculate base pair probabilities for more detailed analysis:

```python
result = vienna.fold('GGGGAAAACCCC', bp_probs=True)
print(f"Number of base pairs: {len(result.bp_probs)}")

# Access base pair probabilities
for bp in result.bp_probs:
    i, j, prob = bp
    if prob > 0.5:  # Only show high-probability pairs
        print(f"Base pair {i}-{j}: probability = {prob:.3f}")
```

### Cofolding Two Sequences

Fold two RNA sequences together to predict their interaction:

```python
# Separate sequences with '&'
result = vienna.cofold('GGGG&AAACCCC')
print(f"Combined structure: {result.dot_bracket}")
print(f"Combined energy: {result.mfe} kcal/mol")
```

### Inverse Folding

Design sequences that fold into a specific structure:

```python
# Target structure and sequence constraints
target_structure = "(((.(((....))).)))"
constraint = "NNNgNNNNNNNNNNaNNN"  # N = any nucleotide, lowercase = specific

# Generate sequences
results = vienna.inverse_fold(target_structure, constraint, n_sol=10)

print(f"Found {len(results)} sequences:")
for seq_score in results:
    print(f"  {seq_score.seq} (score: {seq_score.score:.2f})")
```

### Checking if Sequence Folds to Target Structure

Verify if a sequence folds into a desired structure:

```python
sequence = 'GGGGAAAACCCC'
target = '((((....))))'

if vienna.does_sequence_fold_to(sequence, target):
    print("Sequence folds to target structure!")
else:
    print("Sequence does not fold to target structure")
```

### Working with Longer Sequences

The library handles sequences of any length:

```python
# Longer sequence
long_seq = 'GGGGAAAACCCC' * 10
result = vienna.fold(long_seq)
print(f"Structure length: {len(result.dot_bracket)}")
print(f"Energy: {result.mfe} kcal/mol")
```

### Example: Analyzing a tRNA-like Sequence

```python
# tRNA-like sequence
tRNA_seq = 'GGGGAUAUAGCUCAGUUGGUAGAGCGCUGCCUUUGCACGGCAGAUGUCAGAGGUUCGAUUCUCUGUUAUCCCC'

result = vienna.fold(tRNA_seq, bp_probs=True)
print(f"tRNA structure:\n{result.dot_bracket}")
print(f"Energy: {result.mfe} kcal/mol")
print(f"Ensemble diversity: {result.ens_defect}")

# Find highly probable base pairs
high_prob_pairs = [bp for bp in result.bp_probs if bp[2] > 0.9]
print(f"\nHighly probable base pairs (>90%): {len(high_prob_pairs)}")
```

## Examples

See the `notebooks/` directory for detailed Jupyter notebook examples demonstrating:
- Basic folding workflows
- Base pair probability analysis
- Sequence design with inverse folding
- Cofolding analysis
- And more!