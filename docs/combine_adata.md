### `combine_adata.ipynb` – Detailed Documentation

---

#### 1  Purpose & Position in the Pipeline

`combine_adata.ipynb` **merges all per-study/pre-processed `.h5ad` files into a single AnnData object that is ready for downstream Bayesian modelling**.
Whereas `extract_files.ipynb` produced one `.h5ad` per GEO study (after minimal QC and metadata injection), this notebook:

1. Reads every `*_preprocessed_with_metadata.h5ad` in the working directory
2. Builds the **union set of genes (39 241 features)** across files ([GitHub][1])
3. Re-indexes each AnnData to that union (filling absent genes with 0-sparse columns)
4. Concatenates all objects by rows (cells) with consistent gene order
5. Saves the result as `combined_by_obs.h5ad` (≈ 186 224 cells × 39 241 genes)

This unified object is the direct input for `bayesian_vi.ipynb`, which performs the negative-binomial PyMC model.

---

#### 2  Input Expectations

| Location        | File pattern                        | Origin                                      |
| --------------- | ----------------------------------- | ------------------------------------------- |
| `src/` (or cwd) | `*_preprocessed_with_metadata.h5ad` | Output of `extract_files.ipynb`             |
| Optional        | `gene_lengths.tsv`                  | Already embedded in each `.h5ad` via `.var` |

> **Tip** — Before running, place only the desired studies’ files in the folder or edit the hard-coded `files = [...]` list in **Cell \[#11]**. ([GitHub][1])

---

#### 3  Software Requirements

Same versions as in *extract\_files*, plus **SciPy** for sparse matrices and **h5py** for quick file inspection.

```bash
pip install scanpy anndata scipy h5py
```

---

#### 4  Workflow Breakdown

| Step                                                          | Notebook cell(s) | Details                                                                                                                                                             |
| ------------------------------------------------------------- | ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Imports & file scan**                                       | 1 & 2            | Loads *scanpy*, *anndata*, *glob*, etc.; prints number of `.h5ad` files detected ([GitHub][1])                                                                      |
| **Manual file subset**                                        | 3                | The list `files = [...]` keeps only the six GEO studies used in the project (edit as needed)                                                                        |
| **Gene-union builder**                                        | 4 & 5            | Opens each file in **backed='r'** mode to collect gene indices without loading expression matrices; returns a Python `list` `union_vars` sorted by first-seen order |
| **Re-index helper**                                           | 4                | `reindex_adata()` pads each AnnData with zero-filled sparse columns for missing genes, then re-orders columns to match `union_vars`                                 |
| **Main loop**                                                 | 6                | • Iterates through `files`                                                                                                                                          |
| • Skips corrupt / mislabelled files via `is_valid_h5ad()`     |                  |                                                                                                                                                                     |
| • Applies `reindex_adata()`, then `anndata.concat` (row-wise) |                  |                                                                                                                                                                     |
| • Frees memory with `gc.collect()` after each append          |                  |                                                                                                                                                                     |
| **Save & finish**                                             | final cell       | Writes `combined_by_obs.h5ad` using default compression; prints shape summary to stdout                                                                             |

---

#### 5  Key Implementation Choices

* **Backed mode for union scan** – avoids loading full matrices, cutting peak RAM by >90 %.
* **Sparse zero padding** – retains storage efficiency when adding thousands of absent genes to smaller studies.
* **Duplicate cell IDs** – `.obs_names_make_unique()` warning is printed but not called automatically; run it manually if downstream tools require unique barcodes.
* **Graceful corruption handling** – files lacking an `obs` group are reported and skipped rather than halting execution ([GitHub][1]).

---

#### 6  Running the Notebook

```bash
# activate env
jupyter lab src/combine_adata.ipynb
# Run ▶ “Restart & Run All”
```

Typical runtime: **5–6 min** on a modern laptop (8 GB RAM), peak memory \~4 GB.

---

#### 7  Output Validation

```python
import anndata as ad
adata = ad.read_h5ad("../data/combined_by_obs.h5ad")
print(adata)
# AnnData object with n_obs = 186224 | n_vars = 39241
```

Check that:

* `adata.obs['study']` (or similar) lists exactly six GSE accessions.
* No `nan` in `adata.X.indices` (indicates correct sparse concatenation).

---

#### 8  Troubleshooting & Tips

| Issue                          | Cause                                          | Remedy                                                                                                       |
| ------------------------------ | ---------------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| `Unable to open object 'obs'`  | File is partial or lacks metadata              | Re-run `extract_files.ipynb`, ensuring the conversion step completes                                         |
| MemoryError during concat      | Low RAM + very large combined matrix           | Chunk the concat (e.g., two studies at a time) or use `anndata.read_h5ad(..., backed='r')` with on-disk `.X` |
| Mismatched gene ordering later | Forgot to keep `union_vars` order when slicing | Always slice AnnData columns with `adata = adata[:, union_vars]` after any subsetting                        |

---

#### 9  Next Notebook

Proceed to **`VI_Trace.ipynb`** to:

1. Normalise library sizes (CPM or Scran)
2. Fit hierarchical negative-binomial models in PyMC
3. Infer condition-specific expression changes

[1]: https://raw.githubusercontent.com/HarshiniR4/SCRNAseq_Bayesian_Inference/main/src/combine_adata.ipynb "raw.githubusercontent.com"
