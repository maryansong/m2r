import scanpy as sc
import pandas as pd
import numpy as np

samples = {
    "151510": "151510.h5ad",
    "151673": "151673.h5ad",
    "151676": "151676.h5ad"
}

def preprocess_and_export_pca(filename, sample_id):
    adata = sc.read_h5ad(filename)
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable].copy()
    adata.X = adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    pca_df = pd.DataFrame(adata.obsm['X_pca'], index=adata.obs_names)
    pca_df.to_csv(f"{sample_id}_pca_embeddings.csv")
    print(f"{sample_id} PCA exported.")

for sid, fpath in samples.items():
    preprocess_and_export_pca(fpath, sid)





