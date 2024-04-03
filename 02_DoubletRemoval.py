# Import and Settings

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

# Load Soup-corrected Data

mtx_fltr = !ls ~/AsymDENV/SoupCorrectedGEX | sort  # set working directory as ~/AsymDENV
mtx_fltr = ['SoupCorrectedGEX/' + s for s in mtx_fltr]
mtx_fltr

fltr_ads = [sc.read_10x_mtx(i) for i in mtx_fltr]

sample_id = !ls ~/AsymDENV/SoupCorrectedGEX | sort
sample_id

# Doublet Scoring and Removal

scr_dfs = [sc.external.pp.scrublet(rad, copy=True) for rad in fltr_ads]

scr_dfs[0].obs

for i in range(25):
    fltr_ads[i].obs['scrublet_score'] = scr_dfs[i].obs['doublet_score']
    fltr_ads[i].obs['predicted_doublet'] = scr_dfs[i].obs['predicted_doublet']
fltr_ads[0].obs

raw_ad = anndata.AnnData.concatenate(*fltr_ads, batch_categories=sample_id, index_unique=None)
raw_ad = raw_ad[~raw_ad.obs['predicted_doublet'], :]  # removing doublets

# Write the doublets-removed object
raw_ad.write_h5ad('AsymDENV_raw.h5ad', compression='lzf')