!pip install --quiet git+https://github.com/icbi-lab/scirpy.git@master

%load_ext autoreload
%autoreload 2

import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir
import glob
import os
from matplotlib import pyplot as plt
import seaborn as sns

sc.set_figure_params(figsize=(4, 4))
plt.rcParams['pdf.fonttype'] = 42
sc.settings.verbosity = 2
severity_palette = ['#009966','#ff6600', '#ff0000']
sc.logging.print_header()

# Define the base directory and pattern to match the CSV files (Dandelion output files)
base_dir = "/home/jovyan/AsymDENV/BCR_from_dandelion/"
pattern = "*/dandelion/filtered_contig_igblast_db-pass_genotyped.tsv"

# Use glob to find all files matching the pattern
file_paths = glob.glob(base_dir + pattern)

# Read the files into a list of DataFrames
dfs = []
for file_path in file_paths:
    df = pd.read_csv(file_path, sep='\t')
    sample_id = file_path.split('/')[5]  # Adjust based on your path structure
    df['sequence_id'] = df['sequence_id'].apply(lambda x: sample_id + '_' + x)
    dfs.append(df)

# Concatenate all DataFrames
BCR_combined = pd.concat(dfs, ignore_index=True)
BCR_combined.to_csv("BCR_combined.tsv", sep="\t")


# Load GEX
adata = sc.read_h5ad('AsymDENV_annotated.h5ad')


adata.obs['CellType']


adata = adata[adata.obs['CellType'].isin(['B intermediate', 
                                          'B memory',
                                          'B naive',
                                          'Plasmablast'])]
adata.obs['CellType']


adata.obs['severity'].cat.reorder_categories(['AD', 'DF', 'DHF'], inplace = True)


adata


# Load the merged BCR
adata_BCR = ir.io.read_airr('BCR_merged.tsv')


adata_BCR


ir.pp.merge_with_ir(adata, adata_BCR)


sc.pl.umap(adata, color = 'CellType')


ir.tl.chain_qc(adata)


ax = ir.pl.group_abundance(adata, groupby="receptor_subtype", target_col = 'severity')


ax.figure.set_size_inches(8, 8)
ax.figure.savefig('BCR_QC_IRorNoIR.pdf')


ax = ir.pl.group_abundance(adata, groupby="chain_pairing", target_col = 'severity')


ax.figure.set_size_inches(8, 8)
ax.figure.savefig('BCR_QC_ChainPairing.pdf')


print(
    "Fraction of cells with more than one pair of BCRs: {:.2f}".format(
        np.sum(
            adata.obs["chain_pairing"].isin(
                ["extra VJ", "extra VDJ", "two full chains"]
            )
        )
        / adata.n_obs
    )
)


sc.pl.umap(adata, color="chain_pairing", groups="multichain")


print(
    "Fraction of cells with multichain of BCRs: {:.2f}".format(
        np.sum(
            adata.obs["chain_pairing"].isin(
                ["multichain"]
            )
        )
        / adata.n_obs
    )
)


adata = adata[adata.obs["chain_pairing"] != "multichain", :].copy()


adata = adata[~adata.obs["chain_pairing"].isin(["no IR", "orphan VDJ", "orphan VJ", "ambiguous"]), :].copy()


ax = ir.pl.group_abundance(adata, groupby="chain_pairing", target_col="severity")


# using default parameters, `ir_dist` will compute nucleotide sequence identity
ir.pp.ir_dist(adata)
ir.tl.define_clonotypes(adata, receptor_arms="all", dual_ir="primary_only")


ir.tl.clonotype_network(adata, min_cells=5)


ax = ir.pl.clonotype_network(
    adata, color="donor_id", base_size=30, show_labels=False, panel_size=(9, 9)
)


ax.figure.set_size_inches(12, 10)
ax.figure.savefig('BCR_ClonotypeNW_donorID.pdf')


ax = ir.pl.clonotype_network(
    adata, color="severity", base_size=30, show_labels=False, panel_size=(9, 9), palette = severity_palette
)


ax.figure.set_size_inches(12, 10)
ax.figure.savefig('BCR_ClonotypeNW_severity.pdf')


ir.pp.ir_dist(
    adata,
    metric="alignment",
    sequence="aa",
    cutoff=15,
)


ir.tl.define_clonotype_clusters(
    adata, 
    sequence="aa", 
    metric="alignment", 
   receptor_arms="all", 
   dual_ir="any"
)


ir.tl.clonotype_network(adata, min_cells=10, sequence="aa", metric="alignment")


ax = ir.pl.clonotype_network(
    adata, color="donor_id", show_labels=False, panel_size=(7, 7), base_size=20
)


ax.figure.set_size_inches(12, 10)
ax.figure.savefig('BCR_ClonalGroupNW_donorID.pdf')


ax = ir.pl.clonotype_network(
    adata, color="severity", show_labels=False, panel_size=(7, 7), base_size=20
)


ax.figure.set_size_inches(12, 10)
ax.figure.savefig('BCR_ClonalGroupNW_severity.pdf')


ir.pl.clonotype_network(
    adata, color="s.type", label_fontsize=9, panel_size=(7, 7), base_size=20
)


ir.pl.clonotype_network(
    adata, color="IR_VDJ_1_v_call", show_labels=False, panel_size=(7, 7), base_size=20
)


adata.obs.loc[adata.obs["cc_aa_alignment"] == "737", :].groupby(
    [
        "IR_VJ_1_junction_aa",
        #"IR_VJ_2_junction_aa",
        "IR_VDJ_1_junction_aa",
        #"IR_VDJ_2_junction_aa",
        "receptor_subtype",
    ],
    observed=True,
).size().reset_index(name="n_cells")


group737.to_csv('group737.csv')


ir.tl.clonal_expansion(adata)


sc.pl.umap(adata, color=["clonal_expansion", "clone_id_size"])


sc.set_figure_params(figsize=(8, 8))
ir.pl.clonal_expansion(adata, groupby="CellType", clip_at=4, normalize=False, 
                       fig_kws = {'figsize': (6, 6), 'dpi': 120})


ir.pl.clonal_expansion(adata, "CellType", fig_kws = {'figsize': (6, 6), 'dpi': 120})


adata_pb = adata[adata.obs['CellType'] == 'Plasmablast']


sc.set_figure_params(figsize=(8, 8))
ir.pl.clonal_expansion(adata_pb, groupby="severity", clip_at=4, normalize=True, 
                       fig_kws = {'figsize': (6, 6), 'dpi': 120})


ir.pl.clonal_expansion(adata_pb, "severity", fig_kws = {'figsize': (6, 6), 'dpi': 120})

ax = ir.pl.alpha_diversity(adata, groupby="CellType", fig_kws = {'figsize': (6, 6), 'dpi': 120})

ax = ir.pl.alpha_diversity(adata, groupby="severity", fig_kws = {'figsize': (6, 6), 'dpi': 120})

ax = ir.pl.alpha_diversity(adata_pb, groupby="severity", fig_kws = {'figsize': (6, 6), 'dpi': 120})

ax = ir.pl.alpha_diversity(adata_pb, groupby="donor_id", fig_kws = {'figsize': (6, 6), 'dpi': 120})
#extract data
labels = [l.get_text() for l in ax.get_xticklabels()]
values = [rect.get_height() for rect in ax.patches]
df = pd.DataFrame({'donor_id': labels, 'values': values})

df1 = adata.obs[['donor_id', 'severity']].drop_duplicates(subset=['donor_id'])
df.merge(df1).sort_values(by=['severity']).to_csv('plasmablasts_diversity.csv')

adata_exp = adata[adata.obs['clone_id_size'] > 1]

x = pd.crosstab(adata_exp.obs['CellType'], [adata.obs['severity'], adata.obs['sample_id']], normalize = 'columns')

x.to_csv('prop_clonalB.csv')

x = pd.crosstab(adata.obs['CellType'], [adata.obs['severity'], adata.obs['donor_id']])

x.to_csv('BCR_allB_number.csv')

x = pd.crosstab(adata_exp.obs['CellType'], [adata_exp.obs['severity'], adata_exp.obs['donor_id']])

x.to_csv('BCR_clonalB_number.csv')

ir.pl.group_abundance(adata, groupby="clone_id", target_col="CellType", max_cols=10, fig_kws = {'figsize': (6, 6), 'dpi': 120})

ir.pl.group_abundance(adata, groupby="clone_id", target_col="severity", max_cols=10, fig_kws = {'figsize': (6, 6), 'dpi': 120})

ax = ir.pl.group_abundance(
    adata, groupby="IR_VDJ_1_v_call", target_col="severity", normalize=False,
    max_cols = 10,
    fig_kws = {'figsize': (8, 5), 'dpi': 120}
)

ax.figure.savefig('BCR_Top10_betaV_usage_by_severity_all.pdf')

ax = ir.pl.group_abundance(
    adata_exp, groupby="IR_VDJ_1_v_call", target_col="severity", normalize=False,
    max_cols = 10,
    fig_kws = {'figsize': (8, 5), 'dpi': 120}
)

ax = ir.pl.group_abundance(
    adata_pb, groupby="IR_VDJ_1_v_call", target_col="severity", normalize=False,
    max_cols = 10,
    fig_kws = {'figsize': (8, 5), 'dpi': 120}
)

ax.figure.savefig('BCR_Top10_betaV_usage_by_severity_pb.pdf')

ax = ir.pl.group_abundance(
    adata_pb_exp, groupby="IR_VDJ_1_v_call", target_col="severity", normalize=False,
    max_cols = 10,
    fig_kws = {'figsize': (8, 5), 'dpi': 120}
)

ax = ir.pl.group_abundance(
    adata, groupby="IR_VDJ_1_v_call", target_col="donor_id", normalize=True,
    max_cols = 10,
    fig_kws = {'figsize': (17, 5), 'dpi': 120}
)

ax = ir.pl.vdj_usage(adata, 
                    full_combination=False, max_segments=None, max_ribbons=30,
                    fig_kws = {'figsize': (17, 6), 'dpi': 120})

ax = ir.pl.vdj_usage(adata_pb[adata_pb.obs['severity'] == 'AD'], 
                    full_combination=False, max_segments=None, max_ribbons=30,
                    fig_kws = {'figsize': (17, 6), 'dpi': 120})


adata_exp_pb = adata_exp[adata_exp.obs['CellType'] == 'Plasmablasts']


ax = ir.pl.vdj_usage(adata_exp_pb[adata_exp_pb.obs['severity'] == 'AD'], 
                    full_combination=False, max_segments=None, max_ribbons=30,
                    fig_kws = {'figsize': (17, 6), 'dpi': 120})

df2 = adata_exp_pb.obs[['donor_id', 'severity', 'IR_VJ_1_v_call']]

ax = ir.pl.group_abundance(
    adata, groupby="IR_VDJ_1_c_call", target_col="CellType", normalize=True,
    max_cols = 50,
    fig_kws = {'figsize': (8, 5), 'dpi': 200}
)


ax = ir.pl.group_abundance(
    adata, groupby="IR_VDJ_1_c_call", target_col="severity", normalize=True,
    max_cols = 50,
    fig_kws = {'figsize': (8, 5), 'dpi': 200}
)


ir.pl.vdj_usage(
    adata[adata.obs["cc_aa_alignment"].isin(['805']), :],
    max_ribbons=None,
    max_segments=100,
    fig_kws = {'figsize': (17, 6), 'dpi': 120}
)


ax = ir.pl.vdj_usage(adata, full_combination=False, max_segments=None, max_ribbons=30,
                    fig_kws = {'figsize': (17, 6), 'dpi': 120})


ax = ir.pl.vdj_usage(adata[adata.obs['severity'] == 'AD'], full_combination=False, max_segments=None, max_ribbons=30,
                    fig_kws = {'figsize': (17, 6), 'dpi': 120})


ax = ir.pl.vdj_usage(adata[adata.obs['severity'] == 'DF'], full_combination=False, max_segments=None, max_ribbons=30,
                    fig_kws = {'figsize': (17, 6), 'dpi': 120})


ax = ir.pl.vdj_usage(adata[adata.obs['severity'] == 'DHF'], full_combination=False, max_segments=None, max_ribbons=30,
                    fig_kws = {'figsize': (17, 6), 'dpi': 120})


adata_exp = adata[adata.obs['clone_id_size'] >= 2]
adata_pb = adata[adata.obs['CellType'] == 'Plasmablast']
adata_pb_exp = adata_exp[adata_exp.obs['CellType'] == 'Plasmablast']


test = adata[adata.obs["IR_VDJ_1_v_call"] == 'IGHV3-23*01,IGHV3-23D*01']

ax = ir.pl.group_abundance(
    test, groupby="severity", target_col="CellType")


test = adata[adata.obs["IR_VDJ_1_v_call"] == 'IGHV4-39*01']

ax = ir.pl.group_abundance(
    test, groupby="severity", target_col="CellType")


test = adata[adata.obs["IR_VDJ_1_v_call"] == 'IGHV3-15*01']

ax = ir.pl.group_abundance(
    test, groupby="severity", target_col="CellType")


test = adata_pb[adata_pb.obs["IR_VDJ_1_v_call"] == 'IGHV4-39*01']

ax = ir.pl.group_abundance(
    test, groupby="severity", target_col="CellType")


test = adata_exp[adata_exp.obs["IR_VDJ_1_v_call"] == 'IGHV4-39*01']

ax = ir.pl.group_abundance(
    test, groupby="severity", target_col="CellType")


test = adata_exp[adata_exp.obs["IR_VDJ_2_v_call"] == 'IGHV4-39*01']

ax = ir.pl.group_abundance(
    test, groupby="donor_id", target_col="CellType")

ax = ir.pl.group_abundance(
    adata, groupby='IR_VDJ_1_c_call', target_col="severity", normalize=True,
    max_cols = 20,
    fig_kws = {'figsize': (8, 5), 'dpi': 120}
)


ax.figure.savefig('BCR_Isotype_Normalized.pdf')


#extract data
labels = [l.get_text() for l in ax.get_xticklabels()]
values = [rect.get_height() for rect in ax.patches]
df = pd.DataFrame({'donor_id': labels, 'values': values})

df1 = adata_ig.obs[['donor_id', 'severity']].drop_duplicates(subset=['donor_id'])
df.merge(df1).sort_values(by=['severity']).to_csv('BCR_PB_IgG1_Numbers.csv')
