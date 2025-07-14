
# ITHmapper

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

1. Clone the repository:

   ```bash
   git clone https://github.com/abelson-lab/ITHmapper.git
   cd ITHmapper
   ```

2. Install with pip (recommended for development):

   ```bash
   pip install .
   ```

   (Requires Python 3.7+)

---

## Dependencies

- `louvain`
- `igraph`
- `scanpy`
- `anndata`
- `numpy`
- `pandas`
- `hotspotsc`
- `scikit-learn`
- `scipy`
- `annoy`
- `tqdm`
- `importlib_resources` (for Python <3.9)

All required dependencies are specified in `pyproject.toml` or `setup.py`.

---

## Quick Start

ITHmapper requires an adata file with raw counts in adata.layers["counts"].
Addionally the adata.var.index must be ensembl gene IDs without version (eg. ENSG00000186827, not ENSG00000186827.1).
ITHmapper also calculates HVGs internally, please provide an unfiltered adata file with all genes.
It is okay if HVGs were previously calculated, these will get ignored by the pipeline.
ITHmapper also requires either an scvi or pca embedding of the cells for the Hotspot scoring. 
Ignore the "adata.X seems to be already log-transformed." warning if the input adata was already transformed, ITHmapper is still using the raw counts and re-transforming them, see [scanpy issue](https://github.com/scverse/scanpy/issues/1333).
**Cancer types**
ITHmapper will work with the following cancer_type parameters:

'Bladder', 'Breast', 'Colorectal', 'Gastric',
'Kidney_RCC', 'Liver_CHOL', 'Liver_HCC', 'Lung_LUAD',
'Lung_LUSC', 'Lung_SCC', 'Neuroblastoma', 'Ovarian_HGSOC',
'Pancreas', 'Prostate'

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
    flag_cells = True,
    filter_silhouette = 0.2,
    verbose = True
)

# The mapped consensus state is in:
print(final_adata.obs['cancer_state'])
```

---

## Pipeline Steps

1. **Module Scoring:** Scores query cells for reference Hotspot modules (multiple seeds)
2. **Merging and Scaling:** Merges all module scores and scales/aggregates by meta-module
3. **Neighbor Graph:** Builds a nearest-neighbor graph in meta-module space
4. **Clustering/Filtering:** Clusters cells by Louvain, filters cells by silhouette
5. **Reference Mapping:** Assigns each cell to the closest consensus reference state
6. **Returns:** Filtered, annotated AnnData object

---


## Citing

If you use this pipeline, please cite the relevant preprint or publication (add here).

---

## License

MIT License (see `LICENSE` file)

---

## Contact

For questions or contributions, please contact [Ido Nofech-Mozes](mailto:ido.nofechmozes@mail.utoronto.ca).
