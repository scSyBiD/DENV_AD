!pip install --quiet git+https://github.com/icbi-lab/scirpy.git@master


%load_ext autoreload
%autoreload 2

import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir
from matplotlib import pyplot as plt, cm as mpl_cm
from cycler import cycler


sc.set_figure_params(figsize=(4, 4))
plt.rcParams['pdf.fonttype'] = 42
sc.settings.verbosity = 2


sc.logging.print_header()


# Load GEX
adata = sc.read_h5ad('AsymDENV_annotated.h5ad')


adata = adata[adata.obs['CellType'].isin(['CD4 CTL', 
                                          'CD4 Naive', 
                                          'CD4 Proliferating',
                                          'CD4 TCM',
                                          'CD4 TEM',
                                          'CD8 Naive',
                                          'CD8 Proliferating',
                                          'CD8 TCM',
                                          'CD8 TEM',
                                          'MAIT',
                                          'Treg',
                                          'dnT',
                                          'gdT'
                                         ])]
adata.obs['CellType']


adata.obs['severity'].cat.reorder_categories(['AD', 'DF', 'DHF'], inplace = True)


adata.write_h5ad('AsymDENV_Tcells.h5ad')


# Load TCR
adata_TCR = ir.io.read_10x_vdj('TCR_merged.csv')


adata_TCR


ir.pp.merge_with_ir(adata, adata_TCR)


sc.pl.umap(adata, color = 'CellType')


pd.crosstab(adata.obs['donor_id'], adata.obs['CellType'])


severity_palette = ['#009966','#ff6600', '#ff0000']


ir.tl.chain_qc(adata)


ax = ir.pl.group_abundance(adata, groupby="receptor_subtype", target_col = 'severity')


ax.figure.set_size_inches(6.7, 6.7)
ax.figure.savefig('TCR_QC_IRorNoIR.pdf')


pd.crosstab(adata.obs['severity'], adata.obs["receptor_subtype"])


ax = ir.pl.group_abundance(adata, groupby="chain_pairing", target_col = 'severity')


ax.figure.set_size_inches(8, 8)
ax.figure.savefig('TCR_QC_ChainPairing.pdf')


print(
    "Fraction of cells with more than one pair of TCRs: {:.2f}".format(
        np.sum(
            adata.obs["chain_pairing"].isin(
                ["extra VJ", "extra VDJ", "two full chains"]
            )
        )
        / adata.n_obs
    )
)


sc.pl.umap(adata, color="chain_pairing", groups="multichain")


adata = adata[adata.obs["chain_pairing"] != "multichain", :].copy()


adata = adata[~adata.obs["chain_pairing"].isin(["no IR", "orphan VDJ", "orphan VJ", "ambiguous"]), :].copy()


ax = ir.pl.group_abundance(adata, groupby="chain_pairing", target_col="severity")


pd.crosstab(adata.obs['donor_id'], adata.obs['CellType'])


# using default parameters, `ir_dist` will compute nucleotide sequence identity
ir.pp.ir_dist(adata)
ir.tl.define_clonotypes(adata, receptor_arms="all", dual_ir="primary_only")


ir.tl.clonotype_network(adata, min_cells=10)


ax = ir.pl.clonotype_network(
    adata, color="donor_id", base_size=20, show_labels=False, panel_size=(9, 9)
)


ax.figure.set_size_inches(12, 10)
ax.figure.savefig('TCR_ClonotypeNW_donorID.pdf')


ax = ir.pl.clonotype_network(
    adata, color="severity", base_size=20, show_labels=False, panel_size=(9, 9), palette = severity_palette
)


ax.figure.set_size_inches(12, 10)
ax.figure.savefig('TCR_ClonotypeNW_severity.pdf')


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


ir.tl.clonotype_network(adata, min_cells=20, sequence="aa", metric="alignment")


ir.pl.clonotype_network(
    adata, color="donor_id", label_fontsize=9, panel_size=(7, 7), base_size=20
)


ir.pl.clonotype_network(
    adata, color="severity", label_fontsize=9, panel_size=(7, 7), base_size=20
)


ir.pl.clonotype_network(
    adata, color="s.type", label_fontsize=9, panel_size=(7, 7), base_size=20
)


adata.obs.loc[adata.obs["cc_aa_alignment"] == "76", :].groupby(
    [
        "IR_VJ_1_junction_aa",
        "IR_VJ_2_junction_aa",
        "IR_VDJ_1_junction_aa",
        "IR_VDJ_2_junction_aa",
        "receptor_subtype",
    ],
    observed=True,
).size().reset_index(name="n_cells")


adata.obs.loc[adata.obs["cc_aa_alignment"] == "31", :].groupby(
    [
        "IR_VJ_1_junction_aa",
        "IR_VJ_2_junction_aa",
        "IR_VDJ_1_junction_aa",
        "IR_VDJ_2_junction_aa",
        "receptor_subtype",
    ],
    observed=True,
).size().reset_index(name="n_cells")


ir.tl.clonal_expansion(adata)


ax = sc.pl.umap(adata, color=["clonal_expansion", "clone_id_size"], save = 'TCR_ClonalExpUMAP.pdf')


#sc.set_figure_params(figsize=(8, 8))
ax = ir.pl.clonal_expansion(adata, groupby="CellType", clip_at=4, normalize=False, 
                       fig_kws = {'figsize': (5, 8), 'dpi': 120})
ax.figure.savefig('TCR_ClonalExp_CellType.pdf')


ax = ir.pl.clonal_expansion(adata, "CellType", fig_kws = {'figsize': (6, 6), 'dpi': 120})
ax.figure.savefig('TCR_ClonalExp_CellType_Norm.pdf')


ax = ir.pl.alpha_diversity(adata, groupby="CellType", fig_kws = {'figsize': (6, 6), 'dpi': 120})


ax = ir.pl.alpha_diversity(adata, groupby="severity", fig_kws = {'figsize': (6, 6), 'dpi': 120})


adata_cyto = adata[adata.obs['CellType'] == 'CD4 CTL']
ax = ir.pl.alpha_diversity(adata_cyto, groupby="severity", fig_kws = {'figsize': (6, 6), 'dpi': 120})


ax = ir.pl.alpha_diversity(adata_cyto, groupby="donor_id", fig_kws = {'figsize': (6, 6), 'dpi': 120})
#extract data
labels = [l.get_text() for l in ax.get_xticklabels()]
values = [rect.get_height() for rect in ax.patches]
df = pd.DataFrame({'donor_id': labels, 'values': values})

df1 = adata.obs[['donor_id', 'severity']].drop_duplicates(subset=['donor_id'])
df.merge(df1).sort_values(by=['severity']).to_csv('TCR_CD4_CTL_AlphaDiversity.csv')


adata_exp = adata[adata.obs['clone_id_size'] >= 2]


x = pd.crosstab(adata.obs['CellType'], [adata.obs['severity'], adata.obs['donor_id']])
x


x.to_csv('TCR_T_number.csv')


ir.pl.group_abundance(adata, groupby="clone_id", target_col="CellType", max_cols=10, fig_kws = {'figsize': (6, 6), 'dpi': 120})


ir.pl.group_abundance(adata, groupby="clone_id", target_col="severity", max_cols=10, fig_kws = {'figsize': (6, 6), 'dpi': 120})





adata_MAIT = adata[adata.obs['CellType'] == 'MAIT']


ax = ir.pl.group_abundance(
    adata, groupby="IR_VDJ_1_v_call", target_col="severity",
    max_cols = 10,
    fig_kws = {'figsize': (17, 5), 'dpi': 120}
)


ax.figure.savefig('TCR_Top10_betaV_usage_by_severity.pdf')


ax = ir.pl.group_abundance(
    adata, groupby="IR_VJ_1_j_call", target_col="CellType", normalize=True,
    fig_kws = {'figsize': (17, 5), 'dpi': 120}
)


ax = ir.pl.vdj_usage(adata, full_combination=False, max_segments=None, max_ribbons=30,
                    fig_kws = {'figsize': (17, 6), 'dpi': 120})


ax = ir.pl.vdj_usage(adata[adata.obs['severity'] == 'AD'], full_combination=False, max_segments=None, max_ribbons=30,
                    fig_kws = {'figsize': (17, 6), 'dpi': 120})


ax = ir.pl.vdj_usage(adata[adata.obs['severity'] == 'DF'], full_combination=False, max_segments=None, max_ribbons=30,
                    fig_kws = {'figsize': (17, 6), 'dpi': 120})


ax = ir.pl.vdj_usage(adata[adata.obs['severity'] == 'DHF'], full_combination=False, max_segments=None, max_ribbons=30,
                    fig_kws = {'figsize': (17, 6), 'dpi': 120})


ax = ir.pl.vdj_usage(adata_MAIT, full_combination=False, max_segments=None, max_ribbons=30,
                    fig_kws = {'figsize': (17, 6), 'dpi': 120})
ax.figure.savefig('TCR_MAITCluster_VDJUsage.pdf')


ax = ir.pl.group_abundance(
    adata_MAIT, groupby="severity", target_col="CellType",
    fig_kws = {'figsize': (17, 5), 'dpi': 120}
)


#Filter true MAIT by TRAV1-2
adata_true_MAIT = adata_MAIT[adata_MAIT.obs["IR_VJ_1_v_call"] == 'TRAV1-2']


#Filter true MAIT by TRAJ33/12/20
adata_true_MAIT = adata_true_MAIT[(adata_true_MAIT.obs["IR_VJ_1_j_call"] == 'TRAJ33') | (adata_true_MAIT.obs["IR_VJ_1_j_call"] == 'TRAJ12') | (adata_true_MAIT.obs["IR_VJ_1_j_call"] == 'TRAJ20')]


ax = ir.pl.vdj_usage(adata_true_MAIT, full_combination=False, max_segments=None, max_ribbons=30,
                    fig_kws = {'figsize': (17, 6), 'dpi': 120})
ax.figure.savefig('TCR_TrueMAIT_VDJUsage.pdf')


adata_exp = adata[adata.obs['clone_id_size'] > 2]
adata_exp


#Filter true iNKT by TRAV10 and TRAJ18
adata_true_NKT = adata_MAIT[(adata_MAIT.obs["IR_VJ_1_v_call"] == 'TRAV10') & (adata_MAIT.obs["IR_VJ_1_j_call"] == 'TRAJ18')]


ax = ir.pl.vdj_usage(adata_true_NKT, full_combination=False, max_segments=None, max_ribbons=30,
                    fig_kws = {'figsize': (17, 6), 'dpi': 120})
ax.figure.savefig('TCR_TrueNKT_VDJUsage.pdf')


ax = ir.pl.group_abundance(
    adata_true_NKT, groupby="donor_id", target_col="CellType",
    fig_kws = {'figsize': (17, 5), 'dpi': 120}
)


#Collect data from clonally expanded true MAIT and iNKT
adata_MAIT = adata[adata.obs['CellType'] == 'MAIT']
adata_true_MAIT = adata_MAIT[adata_MAIT.obs["IR_VJ_1_v_call"] == 'TRAV1-2']
adata_true_MAIT = adata_true_MAIT[(adata_true_MAIT.obs["IR_VJ_1_j_call"] == 'TRAJ33') | (adata_true_MAIT.obs["IR_VJ_1_j_call"] == 'TRAJ12') | (adata_true_MAIT.obs["IR_VJ_1_j_call"] == 'TRAJ20')]

#for iNKT
#adata_true_MAIT = adata_MAIT[(adata_MAIT.obs["IR_VJ_1_v_call"] == 'TRAV10') & (adata_MAIT.obs["IR_VJ_1_j_call"] == 'TRAJ18')]

ax = ir.pl.group_abundance(
    adata_true_MAIT, groupby="donor_id", target_col="CellType")

#extract data
labels = [l.get_text() for l in ax.get_xticklabels()]
values = [rect.get_height() for rect in ax.patches]
df = pd.DataFrame({'donor_id': labels, 'values': values})

df1 = adata_true_MAIT.obs[['donor_id', 'severity']].drop_duplicates(subset=['donor_id'])
df.merge(df1).sort_values(by=['severity']).to_csv('TCR_TrueMAIT_all_CellNum.csv')





ir.pl.vdj_usage(
    adata[adata.obs["clone_id"].isin(['1564', '46652', '46678', '46595', '1580', '1593']), :],
    max_ribbons=None,
    max_segments=100,
    fig_kws = {'figsize': (17, 6), 'dpi': 120}
)


ir.pl.vdj_usage(
    adata[adata.obs["clone_id"].isin(['22430', '10271', '33872', '14844']), :],
    max_ribbons=None,
    max_segments=100,
    fig_kws = {'figsize': (17, 6), 'dpi': 120}
)


adata.write_h5ad('AsymDENV_TCR.h5ad')
