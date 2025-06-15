import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# Load data
adata = sc.read_h5ad("151673.h5ad")
clusters = pd.read_csv("151673_cluster_labels.csv", index_col=0)

# Fix the index format 
clusters.index = clusters.index.str.replace(".1", "-1", regex=False)

# Subset to common barcodes to avoid mismatches
common_barcodes = adata.obs_names.intersection(clusters.index)
adata = adata[common_barcodes].copy()
clusters = clusters.loc[common_barcodes]

# Add BayesSpace cluster labels to AnnData
adata.obs["BayesSpace_cluster"] = clusters["cluster"].astype("category")

# Log-transform expression values
sc.pp.log1p(adata)

# Perform differential expression analysis between clusters
sc.tl.rank_genes_groups(
    adata,
    groupby="BayesSpace_cluster",
    method="t-test",
    use_raw=False,
    n_genes=100
)

# Plot top 10 marker genes per cluster and export to PNG
sc.pl.rank_genes_groups(
    adata,
    n_genes=10,
    sharey=False,
    show=False
)
plt.savefig("top_genes_per_cluster.png", dpi=300, bbox_inches="tight")

# Export top 3 marker genes per cluster to CSV
ranked = adata.uns["rank_genes_groups"]
top_genes = {
    group: [ranked["names"][group][i] for i in range(3)]
    for group in ranked["names"].dtype.names
}
pd.DataFrame(top_genes).to_csv("top_3_marker_genes_per_cluster.csv")






