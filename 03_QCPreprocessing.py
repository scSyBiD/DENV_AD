# Single-cell RNAseq of DENV Patients with Varying Severity: QC and Pre-processing
# Preparations

%load_ext autoreload
%autoreload 2

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sn
import anndata
import scanpy as sc

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

np.set_printoptions(linewidth=180)
sc.settings.verbosity = 1

# Processing

raw_ad = sc.read('AsymDENV_raw.h5ad')

# QC and Filtering

rna = raw_ad.copy()

sc.pp.filter_cells(rna, min_genes=200)
sc.pp.filter_genes(rna, min_cells=3)

rna.var['mt'] = rna.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(rna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

rna = rna[rna.obs.pct_counts_mt < 10, :]

sc.pl.violin(rna, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

sc.pl.scatter(rna, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(rna, x='total_counts', y='n_genes_by_counts')

rna = rna[rna.obs.n_genes_by_counts < 5000, :]

rna.write_h5ad('AsymDENV_QC.h5ad')
