import os
import pandas as pd
import scanpy as sc
import functools

# Print Scanpy version and other details
sc.logging.print_header()

# Print the current working directory
print(os.getcwd())

# Load the annotated dataset
adata = sc.read_h5ad('AsymDENV_annotated.h5ad')
print(adata)

# Display the value counts for the 'severity' column
print(adata.obs['severity'].value_counts())

# Replace '/' with '-' in cell type names to avoid path annotation confusion
adata.obs['CellType'] = list(map(lambda st: str.replace(st, "/", "-"), adata.obs['CellType']))
print(adata.obs['CellType'])

# Function to export differential expression results
def exportDEres(adata, key='rank_genes_groups', column=None, filename=None):
    scores = pd.DataFrame(data=adata.uns[key]['scores'][column], index=adata.uns[key]['names'][column])
    lfc = pd.DataFrame(data=adata.uns[key]['logfoldchanges'][column], index=adata.uns[key]['names'][column])
    pvals = pd.DataFrame(data=adata.uns[key]['pvals'][column], index=adata.uns[key]['names'][column])
    padj = pd.DataFrame(data=adata.uns[key]['pvals_adj'][column], index=adata.uns[key]['names'][column])
    try:
        pts = pd.DataFrame(data=adata.uns[key]['pts'][column], index=adata.uns[key]['names'][column])
    except:
        pass
    scores = scores.loc[scores.index.dropna()]
    lfc = lfc.loc[lfc.index.dropna()]
    pvals = pvals.loc[pvals.index.dropna()]
    padj = padj.loc[padj.index.dropna()]
    try:
        pts = pts.loc[pts.index.dropna()]
    except:
        pass
    try:
        dfs = [scores, lfc, pvals, padj, pts]
    except:
        dfs = [scores, lfc, pvals, padj]
    df_final = functools.reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True), dfs)
    try:
        df_final.columns = ['scores', 'logfoldchanges', 'pvals', 'pvals_adj', 'pts']
    except:
        df_final.columns = ['scores', 'logfoldchanges', 'pvals', 'pvals_adj']

    df_final.to_csv(filename, sep='\t')

# Differential expression analysis for each cell type against 'severity', using Wilcoxon test
for x in list(set(adata.obs['CellType'])):
    try:
        print(x)
        adatax = adata[adata.obs['CellType'] == x].copy()
        sc.tl.rank_genes_groups(adatax, groupby='severity', method='wilcoxon', reference='AD', n_genes=40000)
        for y in ['DF', 'DHF']:
            filepath = '08_DGE/wilcoxon'
            if not os.path.exists(filepath):
                os.makedirs(filepath)
            exportDEres(adatax, column=y, filename=filepath + '/' + y + '_' + x + '_wilcoxon.txt')
    except:
        pass

# Differential expression analysis for each cell type against 'severity' with 'DHF' as reference, using Wilcoxon test
for x in list(set(adata.obs['CellType'])):
    try:
        print(x)
        adatax = adata[adata.obs['CellType'] == x].copy()
        sc.tl.rank_genes_groups(adatax, groupby='severity', method='wilcoxon', reference='DHF', n_genes=40000)
        for y in ['AD', 'DF']:
            filepath = '08_DGE/wilcoxon_rev'
            if not os.path.exists(filepath):
                os.makedirs(filepath)
            exportDEres(adatax, column=y, filename=filepath + '/' + y + '_' + x + '_wilcoxon.txt')
    except:
        pass
