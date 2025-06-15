import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad

#Load cluster labels exported from R and raw data
adata = sc.read_h5ad("151673.h5ad")
cluster_df = pd.read_csv("151673_cluster_labels.csv")
cluster_df.set_index("barcode", inplace=True)

# Fix the index format to match adata.obs
# The cluster file uses ".1" while h5ad uses "-1"
cluster_df.index = cluster_df.index.str.replace(".1", "-1", regex=False)

# Merge cluster labels into the AnnData object
adata.obs = adata.obs.merge(cluster_df, left_index=True, right_index=True, how="left")

# Convert to categorical type for plotting
adata.obs["BayesSpace_cluster"] = adata.obs["cluster"].astype("category")

# Check cluster distribution
print("BayesSpace clusters successfully merged:")
print(adata.obs["BayesSpace_cluster"].value_counts())

adata.write("151673_with_BayesSpace_cluster.h5ad")

# Plot settings
genes_to_plot = {
    "ERBB2": "ERBB2_expression.png",
    "IGHG3": "IGHG3_expression.png",
    "CD79A": "CD79A_expression.png",
    "ESR1":  "ESR1_expression.png"
}

vmin_val = 0
vmax_val = 2

# Plot each gene
for gene, filename in genes_to_plot.items():
    if gene in adata.var_names:
        sc.pl.spatial(
            adata,
            color=gene,
            cmap="inferno",
            size=1.5,
            show=False,
            vmin=vmin_val,
            vmax=vmax_val,
        )
        plt.title(f"{gene} Expression", fontsize=12)
        plt.axis("off")
        plt.savefig(filename, dpi=300, bbox_inches="tight")
        plt.close()
    else:
        print(f" {gene} not found in adata.var_names")






