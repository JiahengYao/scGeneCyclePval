import pandas as pd
from cyclum.tuning import CyclumAutoTune
from cyclum.hdfrw import mat2hdf
import matplotlib.pyplot as plt
import os

# Load the expression data
expression_matrix = pd.read_csv('../data/cache/sample_GSM6886536_Sigmoid_03N_CellType_IgA__Plasma_cells/The_compare_shuffle_3/input_data.csv', index_col=0)

# Train the Cyclum model
model = CyclumAutoTune(expression_matrix.values, max_linear_dims=3, epochs=100, rate=2e-4, verbose=20, encoder_width=[30, 20])
model.train(expression_matrix.values, epochs=300, rate=1e-4, verbose=30)
pseudotime = model.predict_pseudotime(expression_matrix.values)
cell_names = expression_matrix.index.tolist()

# Save the predicted pseudotime
pseudotime_df = pd.DataFrame(pseudotime.reshape(-1, 1), index=cell_names, columns=['Pseudotime'])
pseudotime_file = os.path.join('../data/cache/sample_GSM6886536_Sigmoid_03N_CellType_IgA__Plasma_cells/The_compare_shuffle_3', 'pseudotime.h5')
mat2hdf(pseudotime_df, pseudotime_file)

# Retrieve and save the model weights
weights = model.get_weight()
gene_names = expression_matrix.columns.tolist()
feature_names = ['Weights_' + str(i) for i in range(weights.shape[0])]
weights_df = pd.DataFrame(weights, index=feature_names, columns=gene_names)
weights_file = os.path.join('../data/cache/sample_GSM6886536_Sigmoid_03N_CellType_IgA__Plasma_cells/The_compare_shuffle_3', 'weights.h5')
mat2hdf(weights_df, weights_file)

# Plotting is commented out here; uncomment if needed.
# plt.figure(figsize=(10, 5))
# plt.plot(pseudotime)
# plt.title('Pseudotime Analysis')
# plt.savefig(os.path.join('../data/cache/sample_GSM6886536_Sigmoid_03N_CellType_IgA__Plasma_cells/The_compare_shuffle_3', 'pseudotime_plot.png'))
# plt.close()
