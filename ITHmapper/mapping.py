import importlib.resources
import pandas as pd

def load_reference_module(seed: int, cancer_type: str, package_name: str = "ITHmapper") -> pd.Series:
    """
    Load reference Hotspot module assignments packaged with the module.
    
    Returns:
        pd.Series with module membership, indexed by gene symbol.
    """
    filename = f"reference_modules/seed_{seed}early_stopping_fresh_hotspot_scvi_modules_results_{cancer_type}.csv"
    try:
        with importlib.resources.files(package_name).joinpath(filename).open('r') as f:
            df = pd.read_csv(f, index_col=0)
            module_cancer_res = df[['Module']]
            module_cancer_res = module_cancer_res[~module_cancer_res['Module'].isnull()]
            return module_cancer_res['Module']
    except FileNotFoundError:
        raise ValueError(f"Reference module file not found: {filename}")
    



import numpy as np
import scanpy as sc
import hotspot
from anndata import AnnData
from typing import Optional

def score_reference_hotspot_modules(
    adata: AnnData,
    cancer_type: str,
    seed: int = 0,
    embedding_key: str = "X_scVI",
    n_neighbors: int = 30,
    umi_counts_obs_key: str = "nCount_RNA",
    output_csv: Optional[str] = None,
    package_name: str = "ITHmapper"
) -> AnnData:
    """
    Score reference Hotspot modules in query scRNA-seq data.
    Reference modules are loaded automatically from packaged data.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with filtered query cells and embeddings computed, must contain a 'counts' layer with raw counts.
    cancer_type : str
        Cancer type string (e.g. "Lung_LUAD") for selecting reference modules.
    seed : int
        Seed for reproducibility.
    embedding_key : str
        Key in `adata.obsm` for embedding ('X_scVI' or 'X_pca').
    n_neighbors : int
        Number of neighbors for Hotspot KNN graph.
    umi_counts_obs_key : str
        Key in `adata.obs` for total UMI counts.
    output_csv : str, optional
        If given, save output to this CSV.
    package_name : str
        Name of package containing reference module data.
    
    Returns
    -------
    AnnData
        AnnData with module scores in `.obs`
    """
    # Copy to avoid changing original
    adata_rep = adata.copy()
    
    #require that the user provides a layer with counts
    if "counts" not in adata_rep.layers:
        raise ValueError("Input AnnData must contain a 'counts' layer with raw counts data.")
    
    #check if embedding_key is in obsm
    if embedding_key not in adata_rep.obsm:
        raise ValueError(f"Embedding key '{embedding_key}' not found in adata.obsm. Please compute the scvi or pca embedding first.")
    # Normalize, log1p, HVG selection
    adata_rep.X = adata_rep.layers["counts"].copy()
    sc.pp.normalize_total(adata_rep)
    sc.pp.log1p(adata_rep)
    sc.pp.highly_variable_genes(adata_rep, n_top_genes=2000)
    # Load reference modules
    from .mapping import load_reference_module  # If in same module
    module_cancer_res = load_reference_module(seed, cancer_type, package_name=package_name)
    #keep genes that are in module_cancer_res or adata_rep.var.highly_variable
    adata_rep.var['module_gene'] = adata_rep.var.index.isin(module_cancer_res.index)
    adata_rep = adata_rep[:, adata_rep.var['module_gene'] | adata_rep.var['highly_variable']].copy()
    sc.pp.filter_genes(adata_rep, min_cells=3)
    
    # Run Hotspot
    np.random.seed(seed)
    hs = hotspot.Hotspot(
        adata_rep,
        layer_key="counts",
        model='danb',
        latent_obsm_key=embedding_key,
        umi_counts_obs_key=umi_counts_obs_key
    )
    np.random.seed(seed)
    hs.create_knn_graph(weighted_graph=False, n_neighbors=n_neighbors)
    
    # Only keep genes present in both modules and adata
    index_keep = set(module_cancer_res.index) & set(adata_rep.var.index)
    #if no genes match, raise error
    if len(index_keep) == 0:
        raise ValueError(f"No genes from reference module found in query data, are you using the correct format for genes?, used ensembl IDs")
    module_cancer_res = module_cancer_res.loc[list(index_keep)]
    hs.modules = module_cancer_res
    
    # Calculate module scores
    module_scores_cancer_query = hs.calculate_module_scores()
    adata_rep.obs = adata_rep.obs.join(module_scores_cancer_query)
    
    # Optionally save
    if output_csv is not None:
        adata_rep.obs.to_csv(output_csv)
    
    return adata_rep


import importlib.resources
import pandas as pd
import os

def load_consensus_mean_exp_mat(
    cancer_type: str,
    package_name: str = "ITHmapper",
    subdir: str = "reference_modules/consensus_mean_exp_mat",
) -> pd.DataFrame:
    """
    Load consensus mean expression matrix for given cancer type from package data.

    Parameters
    ----------
    cancer_type : str
        The cancer type for which to load the matrix (e.g., 'Lung_LUAD_LUAD').
    package_name : str
        Name of your package.
    subdir : str
        Subdirectory in the package where matrices are stored.

    Returns
    -------
    pd.DataFrame
        DataFrame with gene/module names as rows, reference states as columns.
    """

    
    if cancer_type == "Bladder":
        cancer_type_TCGA = "Bladder_BLCA"
    elif cancer_type == "Brain":
        cancer_type_TCGA = "Brain_GBM"
    elif cancer_type == "Breast":
        cancer_type_TCGA = "Breast_BRCA"
    elif cancer_type == "Colorectal":
        cancer_type_TCGA = "Colorectal_COAD"
    elif cancer_type == "Gastric":
        cancer_type_TCGA = "Gastric_STAD"
    elif cancer_type == "Kidney_RCC":
        cancer_type_TCGA = "Kidney_RCC_KIRC"
    elif cancer_type == "Liver_HCC":
        cancer_type_TCGA = "Liver_HCC_LIHC"
    elif cancer_type == "Lung_LUAD":
        cancer_type_TCGA = "Lung_LUAD_LUAD"
    elif cancer_type == "Neuroblastoma":
        cancer_type_TCGA = "Neuroblastoma_TARGET-NBL"
    elif cancer_type == "Ovarian_HGSOC":
        cancer_type_TCGA = "Ovarian_HGSOC_OV"
    elif cancer_type == "Pancreas":
        cancer_type_TCGA = "Pancreas_PAAD"
    elif cancer_type == "Prostate":
        cancer_type_TCGA = "Prostate_PRAD"

    # Try primary filename
    filename = f"{subdir}/consensus_mean_exp_mat_{cancer_type_TCGA}.tsv"
    with importlib.resources.files(package_name).joinpath(filename).open("r") as f:
        df = pd.read_csv(f, sep=",", index_col=0)
    return df



import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist

def map_cells_to_consensus_states(
    adata,
    meta_scores,
    consensus_mean_exp_mat,
    cancer_type: str,
    flag_cells: bool = True,
    filter_silhouette: float = 0,
):
    """
    Map each filtered cell to the closest consensus state.

    Parameters
    ----------
    adata : AnnData
        Filtered AnnData object (after silhouette filtering).
    meta_scores : pd.DataFrame
        DataFrame of per-cell scaled meta-module scores (rows = cells, cols = meta-modules).
    consensus_mean_exp_mat : pd.DataFrame
        Consensus mean expression matrix (rows = states, cols = meta-modules or vice versa).
    cancer_type : str
        Cancer type string (e.g., "Lung_LUAD") to select correct mapping.
    flag_cells : bool
        If True, flag cells based on silhouette scores before mapping.
        If False, map all cells regardless of silhouette scores.
        Note: This should be set to False if you want to map faster.
    filter_silhouette : float
        Minimum silhouette score to keep a cell.
        If flag_cells is True, cells with silhouette < filter_silhouette will be flagged as unclear.
        If flag_cells is False, this parameter is ignored.
    

    Returns
    -------
    None (modifies adata.obs in place)
    """
    # Align meta_scores to cells in adata
    query_cells = adata.obs_names
    fil_merged_scores = meta_scores.loc[query_cells]
    fil_merged_scores.columns = fil_merged_scores.columns.str.replace('^meta_module_', '', regex=True)
    fil_merged_scores.columns = fil_merged_scores.columns.str.replace('_score', '', regex=True)
    # Ensure columns are in same order as consensus reference
    shared_cols = [c for c in consensus_mean_exp_mat.columns if c in fil_merged_scores.columns]
    fil_merged_scores = fil_merged_scores[shared_cols]
    consensus_mat = consensus_mean_exp_mat[shared_cols]

    # Compute Euclidean distance matrix (cells x reference states)
    dist_centroids = cdist(fil_merged_scores.values, consensus_mat.values, metric="euclidean")

    # Get closest reference state for each cell
    min_distance_state = consensus_mean_exp_mat.index[np.argmin(dist_centroids, axis=1)]

    # Add to AnnData
    adata.obs['min_distance_state'] = min_distance_state
    with importlib.resources.files("ITHmapper").joinpath("reference_modules/data_for_cell_state_dict.csv").open("r") as f:
        mapping_df = pd.read_csv(f)
    # You may need to adjust column names here if not auto-detected:
    mapping_df.columns = ['cancer_type', 'original_label', 'new_label', 'Program Category']
    # Filter to correct cancer type
    mapping_df = mapping_df[mapping_df['cancer_type'] == cancer_type]
    mapping_dict = dict(zip(mapping_df['original_label'], mapping_df['new_label']))
    # Map each cell's label
    adata.obs['cancer_state'] = adata.obs['min_distance_state'].map(mapping_dict).fillna(adata.obs['min_distance_state'])
    mapping_dict_prog_type = dict(zip(mapping_df['original_label'], mapping_df['Program Category']))
    adata.obs['program_type'] = adata.obs['min_distance_state'].map(mapping_dict_prog_type).fillna(adata.obs['min_distance_state'])
    if flag_cells:
        # Flag cells with silhouette < filter_silhouette
        unclear_cells = adata.obs['max_sil_score'] < filter_silhouette
        adata.obs.loc[unclear_cells, 'cancer_state'] = 'unclear'