# Author: jhyao
# Date: Dec 15, 2024
# Description:
# This script processes single-cell gene expression data with pseudotime inference using Cyclum. 
# It generates background distributions by fully shuffling the expression matrix multiple times 
# and subsequently compares the observed results with these null models. This approach facilitates 
# statistical evaluation of gene expression dynamics along pseudotime trajectories.

import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import argparse
import time
import os
import h5py
import matplotlib.pyplot as plt
from cyclum.tuning import CyclumAutoTune
from cyclum.hdfrw import mat2hdf
import itertools

def shuffle_matrix_fully(expression_matrix):
    """ 
    Perform a complete random permutation of all elements in the expression matrix to generate a null model.
    """
    shuffled_array = np.random.permutation(expression_matrix.values.flatten())
    shuffled_matrix = pd.DataFrame(
        shuffled_array.reshape(expression_matrix.shape), 
        index=expression_matrix.index, 
        columns=expression_matrix.columns
    )
    return shuffled_matrix

def check_files_exist(output_base_path, file_prefix):
    """ 
    Check if the pseudotime, weights, and plot files already exist and are non-empty for a given result directory.
    """
    output_dir = os.path.join(output_base_path, file_prefix)
    pseudotime_file = os.path.join(output_dir, 'pseudotime.h5')
    weights_file = os.path.join(output_dir, 'weights.h5')
    pseudotime_plot_file = os.path.join(output_dir, 'pseudotime_plot.png')
    
    return (os.path.exists(pseudotime_file) and os.path.getsize(pseudotime_file) > 0 and
            os.path.exists(weights_file) and os.path.getsize(weights_file) > 0 and
            os.path.exists(pseudotime_plot_file) and os.path.getsize(pseudotime_plot_file) > 0)

def process_and_save_data(expression_matrix, output_dir, file_prefix):
    """ 
    Save input data for cyclum modeling, generate a corresponding Python script for model training, 
    and launch the process in the background. This procedure creates a directory for results and stores input data there.
    """
    output_dir = os.path.join(output_dir, file_prefix)
    os.makedirs(output_dir, exist_ok=True)

    # Save the input expression data for modeling
    input_data_path = os.path.join(output_dir, 'input_data.csv')
    expression_matrix.to_csv(input_data_path)

    # Generate a Python script to train Cyclum, predict pseudotime, and save results
    script_content = f"""import pandas as pd
from cyclum.tuning import CyclumAutoTune
from cyclum.hdfrw import mat2hdf
import matplotlib.pyplot as plt
import os

# Load the expression data
expression_matrix = pd.read_csv('{input_data_path}', index_col=0)

# Train the Cyclum model
model = CyclumAutoTune(expression_matrix.values, max_linear_dims=3, epochs=100, rate=2e-4, verbose=20, encoder_width=[30, 20])
model.train(expression_matrix.values, epochs=300, rate=1e-4, verbose=30)
pseudotime = model.predict_pseudotime(expression_matrix.values)
cell_names = expression_matrix.index.tolist()

# Save the predicted pseudotime
pseudotime_df = pd.DataFrame(pseudotime.reshape(-1, 1), index=cell_names, columns=['Pseudotime'])
pseudotime_file = os.path.join('{output_dir}', 'pseudotime.h5')
mat2hdf(pseudotime_df, pseudotime_file)

# Retrieve and save the model weights
weights = model.get_weight()
gene_names = expression_matrix.columns.tolist()
feature_names = ['Weights_' + str(i) for i in range(weights.shape[0])]
weights_df = pd.DataFrame(weights, index=feature_names, columns=gene_names)
weights_file = os.path.join('{output_dir}', 'weights.h5')
mat2hdf(weights_df, weights_file)

# Plotting is commented out here; uncomment if needed.
# plt.figure(figsize=(10, 5))
# plt.plot(pseudotime)
# plt.title('Pseudotime Analysis')
# plt.savefig(os.path.join('{output_dir}', 'pseudotime_plot.png'))
# plt.close()
"""

    script_path = os.path.join(output_dir, 'run_model.py')
    with open(script_path, 'w') as f:
        f.write(script_content)

    # Launch the script in the background
    os.system(f"nohup python {script_path} > {output_dir}/process.log 2>&1 &")
    print(f"Data processing initiated for {file_prefix}. Check {output_dir}/process.log for progress updates.")
    time.sleep(15)  # Adjust waiting time if required

def process_celltype_data(expression_matrix, metadata, output_prefix, n_shuffles=10, 
                          metadata_filters=None,
                          output_base_path='path_to_save_results'):
    """
    Apply specified filtering conditions to the metadata, select relevant cells, and run Cyclum 
    on both the original and multiple shuffled versions of the expression matrix. 
    
    Parameters:
    expression_matrix (DataFrame): Single-cell gene expression matrix with cells as columns.
    metadata (DataFrame): Corresponding metadata for each cell (index should match cell names).
    output_prefix (str): Prefix for naming output files.
    n_shuffles (int): Number of permutations to generate null distributions.
    metadata_filters (dict): Dictionary specifying filtering criteria {column: [values]}.
    output_base_path (str): Directory path to save the results.
    """
    os.makedirs(output_base_path, exist_ok=True)
    
    # Apply metadata filters
    mask = pd.Series(True, index=metadata.index)
    if metadata_filters:
        for col, values in metadata_filters.items():
            mask &= metadata[col].isin(values)

    selected_cells = metadata[mask].index
    if len(selected_cells) == 0:
        print(f"No cells match the given filters for {output_prefix}: {metadata_filters}. Skipping.")
        return

    selected_expression_matrix = expression_matrix[selected_cells].transpose()

    # Process the original data (if not already done)
    original_file_prefix = output_prefix + "_original"
    if check_files_exist(output_base_path, original_file_prefix):
        print(f"Skipping original data processing as it is already processed for {output_prefix}.")
    else:
        process_and_save_data(selected_expression_matrix, output_base_path, original_file_prefix)

    # Process shuffled data to generate null distributions
    for i in range(n_shuffles):
        file_prefix = output_prefix + f"_shuffle_{i+1}"
        if check_files_exist(output_base_path, file_prefix):
            print(f"Skipping {file_prefix} as it is already processed.")
        else:
            shuffled_matrix = shuffle_matrix_fully(selected_expression_matrix)
            process_and_save_data(shuffled_matrix, output_base_path, file_prefix)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Cyclum pseudotime inference on single-cell data with filtering and null model generation.")
    parser.add_argument('-e', '--expression_matrix', type=str, required=True, help='Path to the input expression matrix CSV file.')
    parser.add_argument('-m', '--metadata', type=str, required=True, help='Path to the input metadata CSV file.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path for output results.')
    parser.add_argument('-p', '--output_prefix', type=str, default='Result', help='Prefix for output directories and files.')
    parser.add_argument('-n', '--n_shuffles', type=int, default=10, help='Number of random permutations to generate null models.')
    # Metadata filters can be specified multiple times, e.g.: -f sample A B -f CellType X Y
    parser.add_argument('-f', '--metadata-filter', nargs='+', action='append', help='Filter metadata: -f ColumnName Value1 Value2 ...')

    args = parser.parse_args()

    # Load input data
    expression_matrix = pd.read_csv(args.expression_matrix, index_col=0)
    metadata = pd.read_csv(args.metadata, index_col=0)

    # Construct metadata_filters dictionary
    metadata_filters = {}
    if args.metadata_filter:
        for filt in args.metadata_filter:
            col = filt[0]
            vals = filt[1:]
            if col not in metadata_filters:
                metadata_filters[col] = vals
            else:
                # Merge and remove duplicates if the same column is specified multiple times
                metadata_filters[col] = list(set(metadata_filters[col] + vals))

    # If multiple filtering criteria are provided with multiple values, generate all combinations
    if metadata_filters:
        keys = list(metadata_filters.keys())
        values_list = [metadata_filters[k] for k in keys]

        for combo in itertools.product(*values_list):
            # Replace spaces and '+' with underscores in directory names
            sanitized_values = [val.replace(' ', '_').replace('+', '_') for val in combo]

            # Create a directory name by joining column-value pairs
            combo_name_parts = []
            for col, val, sval in zip(keys, combo, sanitized_values):
                combo_name_parts.append(f"{col}_{sval}")

            combo_dir_name = "_".join(combo_name_parts)
            output_dir = os.path.join(args.output, combo_dir_name)
            os.makedirs(output_dir, exist_ok=True)

            # Construct a dictionary for the current combination
            combo_filters = {c: [v] for c, v in zip(keys, combo)}

            process_celltype_data(
                expression_matrix=expression_matrix,
                metadata=metadata,
                output_prefix=args.output_prefix,
                n_shuffles=args.n_shuffles,
                metadata_filters=combo_filters,
                output_base_path=output_dir
            )
    else:
        # If no filters are provided, process the entire dataset once
        output_dir = args.output
        os.makedirs(output_dir, exist_ok=True)
        process_celltype_data(
            expression_matrix=expression_matrix,
            metadata=metadata,
            output_prefix=args.output_prefix,
            n_shuffles=args.n_shuffles,
            metadata_filters=None,
            output_base_path=output_dir
        )
