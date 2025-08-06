
# ITHmapper

<p align="center">
    <img src="img/ITHmapper_logo.png" height="240" alt="ITHmapper logo showing stylized interconnected cells forming a network, with vibrant colors suggesting diversity and collaboration, set against a clean white background. The mood is scientific and innovative." style="pointer-events: none;" />
</p>

A Python pipeline for **mapping scRNA-seq query cells to reference cell states** using Hotspot module scoring, meta-module aggregation, consensus state assignment, and cell filtering by clustering quality.


---

## Pipeline Overview

- **Scores query cells** for reference Hotspot modules using a multi-seed approach
- **Merges module scores** and aggregates them into robust meta-modules
- **Scales and summarizes** meta-module scores per cell
- **Builds a neighbor graph** and clusters cells, filtering by silhouette quality
- **Maps each cell** to a consensus reference state by minimum distance in meta-module space
- **Returns** a fully annotated, filtered AnnData object ready for further analysis or visualization

---

## Installation

Install with pip, we strongly recommend to use a new virtual environment:

   ```bash
   pip install --upgrade pip setuptools wheel
   pip install ITHmapper
   ```



## Dependencies

All required dependencies are specified in `pyproject.toml`, to ensure consistent results, versions for all packages are used.

---

## Quick Start

**Input format:**
ITHmapper requires an AnnData object with raw counts in adata.layers["counts"], for large datasets we recommend removing other layers as they will increase the memory used. 
The adata.var.index must be ensembl gene IDs without version (eg. ENSG00000186827, not ENSG00000186827.1).
Please ensure the adata file has all genes and not only HVGs.
ITHmapper also requires either an scvi or pca embedding of the cells for the Hotspot scoring.
ITHmapper requires a key for the adata.obs column containing the number of transcripts in each cell. 
For single dataset/batch samples we recommend using PCA while for samples from many datasets we recommend scVI.
Ignore the "adata.X seems to be already log-transformed." warning if the input adata was already transformed, ITHmapper is still using the raw counts and re-transforming them, see [scanpy issue](https://github.com/scverse/scanpy/issues/1333).

**Cancer types**
Only run ITHmapper on one cancer type at a time.
ITHmapper will work with the following cancer_type parameters, currently other cancer types are not supported:

```python

'Bladder', 'Brain', 'Breast', 'Colorectal', 'Gastric',
'Kidney_RCC', 'Liver_HCC', 'Lung_LUAD',
'Neuroblastoma', 'Ovarian_HGSOC',
'Pancreas', 'Prostate'

```

**Minimal usage:**

```python
import scanpy as sc
from ITHmapper import map_query_to_reference_cell_states

# Load your pre-filtered AnnData (with scVI or PCA embeddings computed)
# the adata must have a 'counts' layer with raw, unnormalized counts.
#the adata.var.index must be ensembl ID eg. ENSG00000186827
adata = sc.read_h5ad("your_filtered_query_cells.h5ad")

# Map to reference cell states
filtered_labelled_adata = map_query_to_reference_cell_states(
    adata,
    cancer_type="Lung_LUAD",
    embedding_key = 'X_scVI',
    umi_counts_obs_key = "nCount_RNA",
    verbose = True
)

# The mapped consensus state is in:
final_adata.obs['cancer_state'].value_counts()
```

---

## Pipeline Parameters

adata: AnnData object to be processed.
cancer_type: one of the cancer types listed above.
embedding_key: embedding for hotspot to use, typically 'X_pca' or 'X_scVI'
flag_cells : bool = False, whether to mark cells that have low silhouette scores as unclear, increases the compute time.
filter_silhouette: float = 0.2, the minimim silhouette score for confident predictions.
umi_counts_obs_key = "nCount_RNA",indicates which column in adata.obs refers to number of transcripts/UMIs.

---


## Citing

If you use this pipeline, please cite the relevant preprint or publication.

---

## License

MIT License (see `LICENSE` file)

---

## Contact

For questions or contributions, please contact [Ido Nofech-Mozes](mailto:ido.nofechmozes@mail.utoronto.ca).
