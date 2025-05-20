import numpy as np
import pandas as pd
import scanpy as sc
import gseapy
import matplotlib.pyplot as plt
import warnings
import os
import subprocess

warnings.simplefilter(action="ignore", category=Warning)

# verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.verbosity = 2

sc.settings.set_figure_params(dpi=80)

path_data = "/your/path/to/data"

path_results = "your/results/"
if not os.path.exists(path_results):
    os.makedirs(path_results, exist_ok=True)

path_file = "./FILE.h5ad"

adata = sc.read_h5ad(path_file)
adata

print(adata.X.shape)
print(type(adata.raw))
print(adata.X[:10,:10])

sc.pl.umap(adata, color='leiden_0.6')

sc.tl.rank_genes_groups(adata, 'leiden_0.6', method='wilcoxon', key_added = "wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key="wilcoxon")

sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, key="wilcoxon", groupby="leiden_0.6", show_gene_labels=True)
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key="wilcoxon", groupby="leiden_0.6")
sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=5, key="wilcoxon", groupby="leiden_0.6")
sc.pl.rank_genes_groups_matrixplot(adata, n_genes=5, key="wilcoxon", groupby="leiden_0.6")

sc.tl.rank_genes_groups(adata, 'leiden_0.6', groups=['1'], reference='2', method='wilcoxon')
sc.pl.rank_genes_groups(adata, groups=['1'], n_genes=20)


sc.pl.rank_genes_groups_violin(adata, groups='1', n_genes=10)

mynames = [x[0] for x in adata.uns['rank_genes_groups']['names'][:10]]
sc.pl.stacked_violin(adata, mynames, groupby = 'leiden_0.6')


cl1 = adata[adata.obs['leiden_0.6'] == '4',:]
cl1.obs['type'].value_counts()

sc.tl.rank_genes_groups(cl1, 'type', method='wilcoxon', key_added = "wilcoxon")
sc.pl.rank_genes_groups(cl1, n_genes=25, sharey=False, key="wilcoxon")

import seaborn as sns

genes1 = sc.get.rank_genes_groups_df(cl1, group='disease', key='wilcoxon')['names'][:5]
genes2 = sc.get.rank_genes_groups_df(cl1, group='Ctrl', key='wilcoxon')['names'][:5]
genes = genes1.tolist() +  genes2.tolist() 
df = sc.get.obs_df(adata, genes + ['leiden_0.6','type'], use_raw=True)
df2 = df.melt(id_vars=["leiden_0.6",'type'], value_vars=genes)

sns.catplot(x = "leiden_0.6", y = "value", hue = "type", kind = 'violin', col = "variable", data = df2, col_wrap=4, inner=None)

genes1 = sc.get.rank_genes_groups_df(cl1, group='Covid', key='wilcoxon')['names'][:5]
genes2 = sc.get.rank_genes_groups_df(cl1, group='Ctrl', key='wilcoxon')['names'][:5]
genes = genes1.tolist() +  genes2.tolist() 

sc.pl.violin(cl1, genes1, groupby='sample', rotation=45)
sc.pl.violin(cl1, genes2, groupby='sample', rotation=45)

genes1 = sc.get.rank_genes_groups_df(cl1, group='disease', key='wilcoxon')['names'][:20]
genes2 = sc.get.rank_genes_groups_df(cl1, group='Ctrl', key='wilcoxon')['names'][:20]
genes = genes1.tolist() +  genes2.tolist() 

sc.pl.dotplot(cl1,genes, groupby='sample')


gene_set_names = gseapy.get_library_name(organism='Human')
print(gene_set_names)

glist = sc.get.rank_genes_groups_df(cl1_sub, group='disease', key='wilcoxon', log2fc_min=0.25, pval_cutoff=0.05)['names'].squeeze().str.strip().tolist()
print(len(glist))

enr_res = gseapy.enrichr(gene_list=glist, organism='Human', gene_sets='GO_Biological_Process_2018', cutoff = 0.5)
enr_res.results.head()


gseapy.barplot(enr_res.res2d,title='GO_Biological_Process_2018')

gene_rank = sc.get.rank_genes_groups_df(cl1_sub, group='disease', key='wilcoxon')[['names','logfoldchanges']]
gene_rank.sort_values(by=['logfoldchanges'], inplace=True, ascending=False)

sc.pp.calculate_qc_metrics(cl1, percent_top=None, log1p=False, inplace=True)

gene_rank = gene_rank[gene_rank['names'].isin(cl1.var_names[cl1.var.n_cells_by_counts>30])]

gene_rank


gene_set_names = gseapy.get_library_name(organism='Human')
print(gene_set_names)

res = gseapy.prerank(rnk=gene_rank, gene_sets='KEGG_2021_Human')

terms = res.res2d.Term
terms[:10]

gseapy.gseaplot(rank_metric=res.ranking, term=terms[0], **res.results[terms[0]])

adata.write_h5ad('./Diff_and_enrich.h5ad')
