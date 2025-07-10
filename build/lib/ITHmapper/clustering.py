import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_samples

def louvain_and_silhouette(
    adata,
    resolutions = np.linspace(0.1, 1.0, 10),
    embedding_key = "X_hotspot",
    neighbors_key = "hotspot_neighbors",
    cluster_prefix = "louvain_",
    silhouette_prefix = "silhouette_",
    filter_silhouette = 0.2,
    random_state = 0,
    verbose = True,
    return_silhouette_matrix = False
):
    """
    Run Louvain clustering at multiple resolutions and compute per-cell silhouette scores.
    Adds max_sil_score to adata.obs and filters cells with max_sil_score < 0.2.
    Optionally returns merged silhouette score DataFrame.

    Parameters
    ----------
    adata : AnnData
        AnnData object with neighbors and embedding computed.
    resolutions : iterable
        List of resolutions to try.
    embedding_key : str
        Key in .obsm for the embedding to compute silhouette scores on.
    neighbors_key : str
        Key prefix in .uns and .obsp for neighbors graph.
    cluster_prefix : str
        Prefix for cluster assignment columns in .obs.
    silhouette_prefix : str
        Prefix for silhouette score columns in .obs.
    filter_silhouette : float
        Minimum silhouette score to keep a cell.
    random_state : int
        For reproducibility.
    verbose : bool
        Print progress.
    return_silhouette_matrix : bool
        If True, also return merged silhouette score DataFrame.

    Returns
    -------
    AnnData (filtered) or (AnnData, silhouette DataFrame)
        Filtered AnnData object, and (optionally) silhouette score matrix.
    """
    sil_matrix = pd.DataFrame(index=adata.obs_names)

    for res in resolutions:
        # Run Louvain clustering
        cluster_key = f"{cluster_prefix}{res:.2f}"
        sc.tl.louvain(
            adata,
            resolution=res,
            key_added=cluster_key,
            random_state=random_state,
            neighbors_key=neighbors_key
        )
        if verbose:
            print(f"Clustering done for resolution {res:.2f}, stored in adata.obs['{cluster_key}']")
        # Compute silhouette scores
        labels = adata.obs[cluster_key].astype(str).values
        embedding = adata.obsm[embedding_key]
        sil_scores = silhouette_samples(embedding, labels, metric='euclidean')
        sil_key = f"{silhouette_prefix}{res:.2f}"
        adata.obs[sil_key] = sil_scores
        sil_matrix[sil_key] = sil_scores
        if verbose:
            print(f"Silhouette scores stored in adata.obs['{sil_key}']")
    
    # Compute and add max_sil_score
    adata.obs['max_sil_score'] = sil_matrix.max(axis=1)
    if verbose:
        print("Appended max_sil_score to adata.obs.")

    # Filter cells with max_sil_score < filter_silhouette
    initial_n = adata.n_obs
    non_confident = adata[adata.obs['max_sil_score'] < filter_silhouette].n_obs

    if verbose:
        print(f"{non_confident} out of {initial_n} cells with max_sil_score < {filter_silhouette}")

    if return_silhouette_matrix:
        # Only return rows for filtered cells
        return adata, sil_matrix.loc[adata.obs_names]
    else:
        return adata