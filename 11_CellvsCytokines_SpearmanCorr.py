import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage

# Load the data
data_path = 'CellDA_cytokine.csv'
data = pd.read_csv(data_path)

# Define cell types and cytokines
cell_types_correct = [
    "ASDC", "B.intermediate", "B.memory", "B.naive", "CD4.CTL", "CD4.Naive",
    "CD4.Proliferating", "CD4.TCM", "CD4.TEM", "CD8.Naive", "CD8.Proliferating",
    "CD8.TCM", "CD8.TEM", "CD14.Mono", "CD16.Mono", "Eryth", "HSPC", "ILC",
    "MAIT", "NK", "NK.Proliferating", "NK_CD56bright", "Plasmablast", "Platelet",
    "Treg", "cDC1", "cDC2", "dnT", "gdT", "pDC"
]

cytokines_correct = [
    "EGF", "Eotaxin", "FGF2", "FLT.3L", "Fractalkine", "G.CSF", "GM.CSF", "GRO",
    "IFN.a2", "IFN.g", "IL.10", "IL.12p40", "IL.12p70", "IL.13", "IL.15",
    "IL.17A", "IL.1alpha", "IL.1beta", "IL.1RA", "IL.2", "IL.4", "IL.5", "IL.6",
    "IL.7", "IL.8", "IP.10", "MCP.1", "MCP.3", "MDC", "MIP.1alpha", "MIP.1beta",
    "PDGF.AA", "TNF.alpha", "TNF.beta", "VEGF"
]

# Filter the data to only include the cell types and cytokines
data_filtered = data[cell_types_correct + cytokines_correct]

# Calculate Spearman's rank correlation coefficients and p-values
spearman_corr_values = pd.DataFrame(index=cell_types_correct, columns=cytokines_correct)
spearman_p_values_matrix = np.zeros((len(cell_types_correct), len(cytokines_correct)))

for i, cell_type in enumerate(cell_types_correct):
    for j, cytokine in enumerate(cytokines_correct):
        corr, p_value = spearmanr(data_filtered[cell_type], data_filtered[cytokine])
        spearman_corr_values.loc[cell_type, cytokine] = corr
        spearman_p_values_matrix[i, j] = p_value

# Convert the correlation coefficients to floats for better visualization
spearman_corr_values = spearman_corr_values.astype(float)

# Calculate the linkage matrices for hierarchical clustering
cell_type_linkage_spearman = linkage(spearman_corr_values, method='ward')
cytokine_linkage_spearman = linkage(spearman_corr_values.transpose(), method='ward')

# Create a mask for p-values < 0.05 to use in the heatmap
significant_mask = spearman_p_values_matrix < 0.05

# Create a clustered heatmap using seaborn's clustermap function with significance annotation
spearman_clustermap_sig = sns.clustermap(spearman_corr_values, 
                                          row_linkage=cell_type_linkage_spearman,
                                          col_linkage=cytokine_linkage_spearman,
                                          figsize=(15, 20),
                                          cmap='coolwarm',
                                          linewidths=.5,
                                          mask=~significant_mask,  # Mask non-significant correlations
                                          dendrogram_ratio=(.1, .2),
                                          cbar_pos=(1, .2, .03, .4))

# Rotate the tick labels for better readability
plt.setp(spearman_clustermap_sig.ax_heatmap.get_xticklabels(), rotation=90)
plt.setp(spearman_clustermap_sig.ax_heatmap.get_yticklabels(), rotation=0)

# Show the plot
plt.show()
