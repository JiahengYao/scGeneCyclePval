# scGeneCyclePval
scGeneCyclePval infers pseudotime from single-cell expression data and evaluates gene expression pattern significance. It uses Cyclum to model trajectories and permutation-based null models from shuffled matrices. Comparing observed and null distributions yields p-values, highlighting significantly dynamic genes.

## Introduction
**scGeneCyclePval** is a computational framework designed for inferring pseudotime dynamics from single-cell expression data and evaluating the significance of observed gene expression patterns. Leveraging the Cyclum model, this tool estimates cellular trajectories and generates background distributions by repeatedly shuffling the original expression matrix to create permutation-based null models. By comparing the empirical distributions from randomized datasets to the observed data, scGeneCyclePval computes p-values and identifies significantly dynamic genes with robust statistical confidence. This approach facilitates a deeper understanding of gene regulatory programs and cellular differentiation processes in complex single-cell transcriptomic datasets.

## Features
- **Pseudotime Inference**: Utilizes the Cyclum model to estimate cellular trajectories.
- **Statistical Validation**: Generates p-values through permutation-based null models to identify significant gene expression patterns.
- **Comparative Analysis**: Supports group-wise comparison of gene weights to identify differences across cell populations.
- **Automated Pipeline**: Provides a shell script to streamline the analysis process.

## Installation

### Prerequisites
- [Anaconda](https://www.anaconda.com/products/individual) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed on your machine.
- Linux operating system.

### 1. Create and Activate Conda Environment
First, create a new Conda environment and activate it:

```bash
conda create -n scGeneCyclePval python=3.7 pip
conda activate scGeneCyclePval
```

### 2. Install Cyclum
Clone the Cyclum repository and install it using pip:
```
git clone https://github.com/KChen-lab/Cyclum.git
cd Cyclum
pip install .
```
During the installation, Cyclum will attempt to install the following dependencies:

python: 3.7.4
keras: 2.2.4
tensorflow: 1.14.0
numpy: 1.16.5
pandas: 0.25.2
scikit-learn: 0.21.3
h5py: 2.9.0
jupyter: 1.0.0
matplotlib: 3.1.1

### 3. Install Additional Python Dependencies
Install the remaining Python libraries required for scGeneCyclePval:
```
conda install pandas numpy scipy h5py matplotlib scikit-learn
```
### 4. Install R Dependencies
The compare.R script requires the argparse package in R. Install it by running the following command in your R environment:
```
install.packages("argparse")
```
5. Clone scGeneCyclePval Repository
Clone the scGeneCyclePval repository to your local machine:
```
git clone https://github.com/JiahengYao/scGeneCyclePval.git
cd scGeneCyclePval
```
**Usage**
Running the Analysis Pipeline
To perform the analysis, execute the provided shell script. This script runs the Python analysis script followed by the R comparative analysis script.
```
bash scripts/ANA_pipline.sh
```
Example Analysis Script
The ANA_pipline.sh script is an example run that does not require any additional parameters. Its content is as follows:
```
#!/bin/bash

# Define input file paths and output directory
expression_matrix="../data/input/expression_matrix.csv"
metadata="../data/input/metadata.csv"
output_dir="../data/cache"

# Define output prefix and number of shuffles
output_prefix="The_compare"
n_shuffles=5

# Run the analysis script
python run_analysis.py \
  --expression_matrix $expression_matrix \
  --metadata $metadata \
  --output $output_dir \
  --output_prefix $output_prefix \
  --n_shuffles $n_shuffles \
  --metadata-filter sample GSM6886536_Sigmoid_03N GSM6886537_Sigmoid_03T \
  --metadata-filter CellType "IgA+ Plasma cells"

# Run the R comparative analysis script
Rscript compare.R --comparison_list ../data/info.txt --output_dir ../data/report1
```
Preparing Info Data
Ensure that your info data is organized as follows:
After running run_analysis.py, prepare an info.txt. This file should contain pairwise comparisons, with each line specifying two groups to compare. The format is as follows:
path_to_group1 path_to_group2

For example:
```
../data/cache/sample_GSM6886536_Sigmoid_03N_CellType_IgA__Plasma_cells ../data/cache/sample_GSM6886537_Sigmoid_03T_CellType_IgA__Plasma_cells
```
Each pair will generate a corresponding subfolder in the --output_dir (e.g., ../data/report1), containing tables of gene cycle effect weight differences, p-values, and FDR values.

Viewing Results
The analysis results will be saved in the specified output directories:

Python Script Outputs (../data/cache):
pseudotime.h5: Pseudotime estimation results.
weights.h5: Gene weights.

R Script Outputs (../data/report1):
Subfolders corresponding to each pairwise comparison listed in info.txt.
Each subfolder contains tables with gene cycle effect weight differences, p-values, and FDR values.

**Contributing**
We welcome contributions! If you have suggestions, bug reports, or want to contribute code, please open an issue or submit a pull request on the GitHub repository.

**License**
This project is licensed under the Apache License 2.0. Additionally, it incorporates Cyclum, which is licensed under the MIT License.

## Acknowledgments
Cyclum: This project utilizes the Cyclum model for pseudotime analysis. We extend its functionality to perform statistical validation of gene expression patterns.
Contributors: We thank all contributors and the open-source community for their valuable tools and resources.

## Dependencies
This tool uses Cyclum for pseudotime analysis, which is licensed under the MIT License.
- Cyclum: https://github.com/KChen-lab/cyclum



