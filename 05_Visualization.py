# Import and Settings

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import anndata
import scanpy as sc
import re
from adjustText import adjust_text

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

np.set_printoptions(linewidth=180)
sc.settings.verbosity = 1

!pip install --quiet harmonypy

# Import the post QC h5ad file
adata = sc.read('AsymDENV_QC.h5ad')

# Set the random seed for reproducibility
np.random.seed(42)

# Normalization (in case some downstream approaches required log-normalized counts)
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Find variable genes (just in case)
sc.pp.highly_variable_genes(adata, n_top_genes=3000, layer="counts", flavor="seurat_v3", batch_key="batch")
sc.pl.highly_variable_genes(adata)

# Remove BCR and TCR genes from highly variable genes
for i in adata.var.index:
    if re.search('^IG[HKL][VDJC]', i) or re.search('^TR[AB][VDJC]', i):
        adata.var.at[i, 'highly_variable'] = False
    if re.search('TRAV1-2', i):
        adata.var.at[i, 'highly_variable'] = True
sc.pl.highly_variable_genes(adata)

# Add metadata
meta = pd.read_csv('meta_full.csv')
adata.obs = adata.obs.reset_index().merge(meta, how="left").set_index('index')

# Add MMRM celltype annotation
celltype = pd.read_csv('MMRM_celltype.csv', index_col=0)
adata.obs = adata.obs.join(celltype)

# Add spca and umap
spca = pd.read_csv('MMRM_spca.csv', index_col=0)
sumap = pd.read_csv('MMRM_umap.csv', index_col=0)
adata.obsm['X_pca'] = spca.to_numpy()
adata.obsm['X_umap'] = sumap.to_numpy()

# Define a palette for the plots
palette = ["#9ca557", "#9f50c8", "#83bd3b", "#6b6bdb", "#beb134", "#df6ed0", "#49c269", "#d33990", "#4e8f2d", 
           "#9b4e9a", "#76b876", "#d54668", "#51c2ac", "#d5433c", "#4cc3e3", "#c55a26", "#628bd8", "#dc9230",
           "#655ea3", "#cca55e", "#bb91da", "#3b7e44", "#df85ab", "#30866c", "#9c446a", "#606c23", "#4b97c8",
           "#e38e6c", "#90682d", "#af5954"]

# Generate labels for UMAP plots
def gen_mpl_labels(adata, groupby, exclude=(), ax=None, adjust_kwargs=None, text_kwargs=None):
    if adjust_kwargs is None:
        adjust_kwargs = {"arrowprops": {"arrowstyle": "-", "color": "k"}, "ax": ax}
    if text_kwargs is None:
        text_kwargs = {"weight": "bold", "ha": "center", "va": "center", "fontsize": 8}

    medians = {}
    for g, g_idx in adata.obs.groupby(groupby).groups.items():
        if g in exclude:
            continue
        medians[g] = np.median(adata.obsm["X_umap"][g_idx], axis=0)

    texts = [ax.text(medians[g][0], medians[g][1], g, **text_kwargs) for g in medians if g not in exclude]
    adjust_text(texts, **adjust_kwargs)

# Subsampling for 'donor_id', 'sex', and 'age'
# Define number of cells to sample from each group
n_cells_per_batch = 10000  # Adjust this to your desired number or to the smallest group size

# Function to subsample and plot UMAP for a given category
def subsample_and_plot(adata, category, n_cells_per_batch, palette, save_filename):
    sampled_indices = pd.DataFrame()

    groups = adata.obs[category].unique()
    for group in groups:
        group_indices = adata.obs.index[adata.obs[category] == group]
        n_samples = min(n_cells_per_batch, len(group_indices))
        sampled_group_indices = np.random.choice(group_indices, size=n_samples, replace=False)
        sampled_indices = sampled_indices.append(pd.DataFrame(sampled_group_indices, columns=['index']))

    adata_sampled = adata[sampled_indices['index'].values, :]
    sc.pl.umap(adata_sampled, color=category, palette=palette, size=10, alpha=0.7, edgecolor='none', show=True, save=save_filename)

# Subsample and plot UMAP for 'donor_id', 'sex', and 'age'
subsample_and_plot(adata, 'donor_id', n_cells_per_batch, palette, 'UMAP_donors.pdf')
subsample_and_plot(adata, 'sex', n_cells_per_batch, palette, 'UMAP_sex.pdf')
subsample_and_plot(adata, 'age', n_cells_per_batch, palette, 'UMAP_age.pdf')

# Cell type proportion analysis
def get_cluster_proportions(adata, cluster_key="cluster_final", sample_key="replicate"):
    sizes = adata.obs.groupby([cluster_key, sample_key]).size()
    props = sizes.groupby(level=1).apply(lambda x: 100 * x / float(x.sum()))
    return props.unstack()

props = get_cluster_proportions(adata, cluster_key="CellType", sample_key="donor_id")
props.plot(kind='bar', stacked=True, legend=False, figsize=(10, 6))
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel('% of cells')
plt.show()

# Write the annotated object to a file
adata.write_h5ad('AsymDENV_annotated.h5ad')
