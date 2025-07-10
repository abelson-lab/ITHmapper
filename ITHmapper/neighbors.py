import scanpy as sc
import numpy as np

def add_hotspot_neighbors(
    adata,
    meta_scores,
    obsm_key="X_hotspot",
    k_param=20,
    metric="euclidean",
    verbose=True
):
    """
    Appends scaled meta-module scores to adata.obsm[obsm_key] and computes a neighbor graph.

    Parameters
    ----------
    adata : AnnData
        Your AnnData object.
    meta_scores : pd.DataFrame
        DataFrame (cells x meta-modules) of scaled mean scores, index matches adata.obs_names.
    obsm_key : str
        Key to store in adata.obsm.
    k_param : int
        Number of neighbors (like Seurat's k.param).
    metric : str
        Distance metric (Seurat default is 'euclidean').
    verbose : bool
        Print progress.

    Returns
    -------
    None (updates adata in place).
    """
    if verbose:
        print(f"Adding {meta_scores.shape[1]}-dimensional X_hotspot to adata.obsm['{obsm_key}']")
    # Align index
    meta_scores = meta_scores.loc[adata.obs_names]
    adata.obsm[obsm_key] = meta_scores.values.astype(np.float32)

    # Set up scanpy neighbors parameters
    neighbors_kwargs = dict(
        n_neighbors=k_param,
        use_rep=obsm_key,
        metric=metric,
        key_added=f"hotspot_neighbors",
    )
    # Compute neighbors
    sc.pp.neighbors(adata, **neighbors_kwargs)
    if verbose:
        print(f"Neighbor graph stored in adata.obsp['hotspot_neighbors_distances'] and ['hotspot_neighbors_connectivities']")
