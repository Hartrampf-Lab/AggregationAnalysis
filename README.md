# Aggregation Analysis
Aggregation Analysis is a tool for investigating and detecting aggregation for the Automated Fast-Flow Peptide Synthesiser. 
It contains the code to implement the methods and concepts introduced in our paper:
<br>
<br>
**A robust analytical method to investigate sequence dependence in flow-based peptide synthesis** 
<br>
(doi: tba) 
<br>
## Conda environment
A conda environment YML file is included in order to successfully run the code.
Create the dataanlysis environment by running the following command: 
<br>
```bash
conda env create -f environment.yml
```
In order to be able to use the environement kernel on jupyter notebook, run the following command:
```bash
conda activate dataanalysis
python -m ipykernel install --user --name dataanalysis
```

## Usage

## Requirements
* <a href='https://www.python.org/downloads/release/python-3110/'>Python 3.11</a>
* <a href='https://www.rdkit.org/'>RDKit </a>
...
