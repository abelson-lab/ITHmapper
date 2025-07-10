
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
    filename = f"reference_modules/cell_state_dictionary_{cancer_type}.json"
    with importlib.resources.files(package_name).joinpath(filename).open("r") as f:
        cell_state_dict = json.load(f)

    # Compute mean of z-scored module columns per meta-module
    meta_scores = pd.DataFrame(index=merged_scores_per_cell.index)
    for meta_idx, module_cols in cell_state_dict.items():
        # Only use columns that exist in the DataFrame
        cols_to_use = [col for col in module_cols if col in merged_scores_per_cell.columns]
        if cols_to_use:
            # Z-score scale columns (axis=0: column-wise)
            scaled = merged_scores_per_cell[cols_to_use].apply(lambda x: (x - x.mean()) / x.std(), axis=0)
            meta_scores[f"{output_prefix}{meta_idx}_score"] = scaled.mean(axis=1)
        else:
            meta_scores[f"{output_prefix}{meta_idx}_score"] = float("nan")
    return meta_scores