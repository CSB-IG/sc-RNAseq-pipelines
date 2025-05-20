# Core scverse libraries
import scanpy as sc
import anndata as ad
sc.settings.set_figure_params(dpi=50, facecolor="white")


path_results = "your/results/"
if not os.path.exists(path_results):
    os.makedirs(path_results, exist_ok=True)

path_file = "./preprocess_FILE.h5ad"
if fetch_data and not os.path.exists(path_file):
    file_url = os.path.join(path_data, "./FILE.h5ad")
    subprocess.call(["curl", "-u", curl_upass, "-o", path_file, file_url ])

adata = sc.read_h5ad(path_file)

for res in [0.02, 0.5, 2.0]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    )

sc.pl.umap(
    adata,
    color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
    legend_loc="on data",
)


sc.pl.umap(
    adata,
    color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
    legend_loc="on data",
)

marker_genes = { HERO GOES THE LIST OF YOUR MARKER GENES

EXAMPLE:

"CD14+ Mono": ["FCN1", "CD14"],
"CD16+ Mono": ["TCF7L2", "FCGR3A", "LYN"],

}

sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.02", standard_scale="var")

adata.obs["cell_type_lvl1"] = adata.obs["leiden_res_0.02"].map(
    {
        "0": "Lymphocytes", #THIS IS AN EXAMPLE OF HOW TO NAME CLUSTERS
        "1": "Monocytes",
        "2": "Erythroid",
        "3": "B Cells",
    }
)

sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.50", standard_scale="var")

FAST DIFFERENTIAL EXPRESSION ANALYSIS

# Obtain cluster-specific differentially expressed genes
sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.50", method="wilcoxon")

sc.pl.rank_genes_groups_dotplot(
    adata, groupby="leiden_res_0.50", standard_scale="var", n_genes=5
)


sc.get.rank_genes_groups_df(adata, group="7").head(5)

dc_cluster_genes = sc.get.rank_genes_groups_df(adata, group="7").head(5)["names"]
sc.pl.umap(
    adata,
    color=[*dc_cluster_genes, "leiden_res_0.50"],
    legend_loc="on data",
    frameon=False,
    ncols=3,
)


adata.write_h5ad("./annot_FILE.h5ad")
