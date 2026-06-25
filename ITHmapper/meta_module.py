
import pandas as pd
import importlib.resources
import json

def compute_scaled_meta_module_scores_for_cancer_type(
    merged_scores_per_cell: pd.DataFrame,
    cancer_type: str,
    package_name: str = "ITHmapper",
    output_prefix: str = "meta_module_"
) -> pd.DataFrame:
    """
    Compute mean meta-module scores (scaled/z-scored per module) for each cell for the specified cancer type.

    Parameters
    ----------
    merged_scores_per_cell : pd.DataFrame
        DataFrame of merged module scores with one row per cell and module score columns.
    cancer_type : str
        Cancer type string (e.g., "Lung_LUAD") to select cell state dictionary.
    package_name : str
        Name of package containing the reference_modules data directory.
    output_prefix : str
        Prefix for output meta-module columns.

    Returns
    -------
    pd.DataFrame
        DataFrame with one column per meta-module, indexed by cell, with mean *scaled* scores.
    """
   # Load cell state dictionary from package data
    dict_filename = f"reference_modules/cell_state_dictionary_{cancer_type}.json"
    with importlib.resources.files(package_name).joinpath(dict_filename).open("r") as f:
        cell_state_dict = json.load(f)

    # 2. Load pre-calculated reference scaling parameters (means and stds)
    scaling_filename = f"reference_modules/cell_state_scaling_{cancer_type}.json"
    with importlib.resources.files(package_name).joinpath(scaling_filename).open("r") as f:
        reference_scaling = json.load(f)

    # Compute mean of reference-scaled module columns per meta-module
    meta_scores = pd.DataFrame(index=merged_scores_per_cell.index)
    for meta_idx, module_cols in cell_state_dict.items():
        # Only use columns that exist in the DataFrame
        cols_to_use = [col for col in module_cols if col in merged_scores_per_cell.columns]
        if cols_to_use:
            # Scaled using reference parameters instead of query statistics
            scaled = merged_scores_per_cell[cols_to_use].copy()
            for col in cols_to_use:
                if col in reference_scaling:
                    ref_mean = reference_scaling[col]["mean"]
                    ref_std = reference_scaling[col]["std"]
                else:
                    # Fallback to query statistics if missing
                    ref_mean = merged_scores_per_cell[col].mean()
                    ref_std = merged_scores_per_cell[col].std()
                
                # Prevent division by zero if std is somehow 0
                ref_std = max(ref_std, 1e-6)
                scaled[col] = (scaled[col] - ref_mean) / ref_std
            
            meta_scores[f"{output_prefix}{meta_idx}_score"] = scaled.mean(axis=1)
        else:
            meta_scores[f"{output_prefix}{meta_idx}_score"] = float("nan")
    return meta_scores