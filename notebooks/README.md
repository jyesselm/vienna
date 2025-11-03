# Vienna Notebook Examples

This directory contains Jupyter notebook examples demonstrating various features of the Vienna package.

## Notebooks

1. **01_basic_folding.ipynb** - Introduction to basic RNA folding
   - Simple sequence folding
   - Visualizing structures
   - Comparing multiple sequences
   - Energy analysis

2. **02_base_pair_probabilities.ipynb** - Base pair probability analysis
   - Calculating base pair probabilities
   - Visualizing probability matrices
   - Analyzing high-probability pairs
   - Probability distributions

3. **03_cofolding.ipynb** - RNA-RNA interaction prediction
   - Basic cofolding
   - Comparing individual vs cofolded structures
   - Interaction energy analysis
   - Testing different sequence combinations

4. **04_inverse_folding.ipynb** - Sequence design
   - Designing sequences for target structures
   - Working with constraints
   - Verifying designs
   - Score analysis

## Running the Notebooks

To run these notebooks, you'll need:

1. Install the Vienna package (see main README)
2. Install Jupyter:
   ```bash
   pip install jupyter matplotlib numpy
   ```
3. Start Jupyter:
   ```bash
   jupyter notebook
   ```
4. Navigate to the `notebooks/` directory and open any notebook

## Requirements

The notebooks require:
- `vienna` package (installed)
- `matplotlib` for plotting
- `numpy` for numerical operations
- `jupyter` to run the notebooks

All dependencies can be installed via:
```bash
pip install jupyter matplotlib numpy
```

Or if using conda:
```bash
conda install jupyter matplotlib numpy
```

