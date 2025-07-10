import pandas as pd
from .mapping import score_reference_hotspot_modules

def merge_hotspot_module_scores_over_seeds(
    adata,
    cancer_type,
    seeds_use=range(1, 11),
    embedding_key="X_scVI",
    umi_counts_obs_key = "nCount_RNA",
    output_csv=None
):
    merged_scores_per_cell = None

    for seed in seeds_use:
        print(f"Scoring modules for seed {seed}...")
        adata_scored = score_reference_hotspot_modules(
            adata,
            cancer_type=cancer_type,
            seed=seed,
            embedding_key=embedding_key,
            umi_counts_obs_key=umi_counts_obs_key,
            output_csv=None  # don't save each, we'll merge first
        )
        # Extract only module score columns: assume these are not present in original obs
        # Exclude any metadata columns you don't want (e.g., obs columns from input)
        adata_scored.obs.columns = adata_scored.obs.columns.astype(str)
        module_cols = [col for col in adata_scored.obs.columns if col.endswith('.0')]
        scores_per_cell = adata_scored.obs[module_cols].copy()
        # Rename columns to include seed
        scores_per_cell.columns = [
            f"seed___{seed}___cancer_cell_module_{col}".replace(".0", "") for col in scores_per_cell.columns
        ]
        # Add cell index if not present
        scores_per_cell.index.name = 'cell'
        scores_per_cell.reset_index(inplace=True)
        
        if merged_scores_per_cell is None:
            merged_scores_per_cell = scores_per_cell
        else:
            # Merge on cell index
            merged_scores_per_cell = pd.merge(
                merged_scores_per_cell, scores_per_cell, on="cell", how="outer"
            )

    if output_csv is not None:
        merged_scores_per_cell.to_csv(output_csv, index=False)

    return merged_scores_per_cell
