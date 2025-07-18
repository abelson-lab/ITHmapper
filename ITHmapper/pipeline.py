import scanpy as sc
from .scoring import merge_hotspot_module_scores_over_seeds
from .meta_module import compute_scaled_meta_module_scores_for_cancer_type
from .neighbors import add_hotspot_neighbors
from .clustering import louvain_and_silhouette
from .mapping import load_consensus_mean_exp_mat, map_cells_to_consensus_states
import numpy as np
from anndata import AnnData

def map_query_to_reference_cell_states(
    adata: AnnData,
    cancer_type: str,
    package_name: str = "ITHmapper",
    embedding_key: str = "X_scVI",
    k_param: int = 20,
    seeds_use = range(1, 11),
    flag_cells : bool = False,
    filter_silhouette: float = 0.2,
    resolutions = np.linspace(0.1, 1.0, 10),
    umi_counts_obs_key = "nCount_RNA",
    verbose: bool = True
) -> AnnData:
    """
    Complete one-liner pipeline to map query cells to reference cell states.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object (filtered to query cells, embedding computed).
    cancer_type : str
        Cancer type string for mapping, must be one of:
        'Bladder', 'Brain','Breast', 'Colorectal', 'Gastric',
        'Kidney_RCC', 'Liver_HCC', 'Lung_LUAD',
        'Neuroblastoma', 'Ovarian_HGSOC',
        'Pancreas', 'Prostate'.
    package_name : str
        Package name for meta-modules and module reference data.
    embedding_key : str
        Key in .obsm for the embedding to use.
    k_param : int
        k for nearest neighbors.
    seeds_use : iterable
        Seeds for Hotspot runs.
    flag_cells : bool
        Whether to flag cells by silhouette score.
        If True, cells with silhouette < filter_silhouette will be flagged as unclear.
    filter_silhouette : float
        Minimum silhouette for filtering.
    resolutions : iterable
        Resolutions for Louvain clustering.
    umi_counts_obs_key : str
        Key in .obs for UMI counts, required for hotspot (e.g., "nCount_RNA").
    verbose : bool
        Verbosity.

    Returns
    -------
    AnnData
        Final filtered and annotated AnnData object (with .obs['cancer_state']).
    """

    ##cancer_type can be:
    # 'Bladder', 'Brain', 'Breast', 'Colorectal', 'Gastric', 'Kidney', 'Liver_HCC','Lung_LUAD',
    # 'Neuroblastoma', 'Ovarian_HGSOC', 'Pancreas', 'Prostate'
    # check and provide error if cancer_type is not in the list
    valid_cancer_types = [
        'Bladder', 'Brain','Breast', 'Colorectal', 'Gastric',
        'Kidney_RCC','Liver_HCC', 'Lung_LUAD',
        'Neuroblastoma', 'Ovarian_HGSOC',
        'Pancreas', 'Prostate'
    ]
    if cancer_type not in valid_cancer_types:
        raise ValueError(
            f"Invalid cancer_type '{cancer_type}'. "
            f"Must be one of: {', '.join(valid_cancer_types)}"
        )

    # 1. Score reference Hotspot modules for all seeds and merge
    merged_scores = merge_hotspot_module_scores_over_seeds(
        adata,
        cancer_type,
        seeds_use=seeds_use,
        embedding_key=embedding_key,
        umi_counts_obs_key=umi_counts_obs_key,
        output_csv=None
    )

    # 2. Compute mean (scaled) meta-module scores
    meta_scores = compute_scaled_meta_module_scores_for_cancer_type(
        merged_scores,
        cancer_type=cancer_type,
        package_name=package_name
    )

    meta_scores.index = merged_scores['cell']

    if flag_cells:
        # 3. Compute initial nearest neighbors using meta-module scores
        add_hotspot_neighbors(
            adata,
            meta_scores,
            obsm_key="X_hotspot",
            k_param=k_param,
            metric="euclidean",
            verbose=verbose
            )
        # 4. Louvain & silhouette, flag by silhouette
        adata, silhouette_score_merged = louvain_and_silhouette(
            adata,
            resolutions=resolutions,
            embedding_key="X_hotspot",
            neighbors_key="hotspot_neighbors",
            return_silhouette_matrix=True,
            filter_silhouette=filter_silhouette,
            verbose=verbose
        )

    # 5. Load consensus mean expression matrix for mapping
    ref_mat = load_consensus_mean_exp_mat(
        cancer_type=cancer_type,
        package_name=package_name,
    )

    # 6. Map filtered cells to reference consensus states
    map_cells_to_consensus_states(
        adata,
        meta_scores,
        ref_mat, 
        flag_cells=flag_cells,
        filter_silhouette=filter_silhouette,
        cancer_type=cancer_type
    )

    if verbose:
        print("Pipeline complete. Return filtered AnnData with mapped consensus cell state labels.")

    return adata

