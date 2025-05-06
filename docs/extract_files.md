### `extract_files.ipynb` – Detailed Documentation

*(part of the **SCRNAseq\_Bayesian\_Inference** project)*

---

#### 1  Purpose & Context

Single-cell RNA-seq (scRNA-seq) data for microglia in Alzheimer’s disease are scattered across several publicly available GEO studies.
`extract_files.ipynb` automates **one-time data extraction and harmonisation** of six GEO Series (GSE) accessions:

| Accession | Tissue / Model                     | Cells   | Notes                       |
| --------- | ---------------------------------- | ------- | --------------------------- |
| GSE98969  | AD mouse cortex (10X)              | 9 ,197  | 2017 Keren-Shaul *et al.*   |
| GSE98971  | AD mouse hippocampus (10X)         | 8 ,638  | Paired with GSE98969        |
| GSE121654 | Human frontal cortex (Drop-seq)    | 70 ,634 | Pre-processed matrix public |
| GSE123025 | Human parietal cortex (10X)        | 1 ,922  | Myeloid-enriched            |
| GSE134379 | Human temporal cortex (Smart-seq2) | 2 ,738  | Full-length protocol        |
| GSE150430 | AD mouse cerebellum (10X)          | 3 ,412  | Late-stage model            |

The notebook yields a **single, concatenated `AnnData` object** (`combined_by_obs.h5ad`) that feeds directly into the downstream Bayesian inference workflow implemented in `src/bayesian_vi.ipynb`.

---

#### 2  Input Files & Directory Structure

```text
data/
├── raw/
│   ├── GSE98969_RAW.tar
│   ├── GSE98971_RAW.tar
│   ├── ...
│   └── GSE150430_RAW.tar
├── meta/
│   ├── GSE98969_SraRunTable.csv
│   ├── ...
│   └── gene_lengths.tsv          # Ensembl gene length table
└── extracted/                    # Populated by the notebook
```

* **`*_RAW.tar`** – the GEO “Supplementary files” archives (downloaded manually or via `wget`).
* **`SraRunTable.csv` / `metadata.tsv`** – GEO’s sample metadata; column names vary by study.
* **`gene_lengths.tsv`** – lengths (bp) for NB regression offsets in later modelling.

Place all archives under `data/raw/` before running the notebook. No files are committed to the repo; `.gitignore` excludes all `*.tar` and `*.h5ad` artefacts.

---

#### 3  Software Requirements

| Package                               | Min. version | Purpose                   |
| ------------------------------------- | ------------ | ------------------------- |
| **Python ≥3.10**                      | –            | interpreter               |
| `scanpy`                              | 1.10         | AnnData handling          |
| `anndata`                             | 0.10         | I/O layer                 |
| `pandas`, `numpy`                     | latest       | tabular ops               |
| `tarfile`, `gzip`, `shutil` (std-lib) | –            | archive extraction        |
| `louvain`                             | 0.8          | invoked for QC clustering |

Install with:

```bash
pip install scanpy anndata pandas louvain
```

*(Conda-based environments are recommended to pin HDF5 and system-level deps.)*

---

#### 4  Workflow Overview

```
   ┌──Download GEO archives
   │
   ├─▶ 1. Untar each *_RAW.tar → sample-level .gz files
   │
   ├─▶ 2. Decompress *.gz → plain-text count matrices (*.txt)
   │
   ├─▶ 3. Parse GEO metadata tables  ┐
   │                                 │
   ├─▶ 4. Load counts into pandas    ├──▶ 5. Harmonise gene IDs & sample
   │                                 │       barcodes across studies
   └─▶ 6. Build AnnData per study    ┘
                  │
                  ▼
   7. Concatenate AnnData objects (outer join)
   8. Inject unified metadata (study, sample, organism, technology, condition)
   9. Quality-control sanity checks (cell/gene counts, mitochondrial %, etc.)
  10. Save as combined_by_obs.h5ad
```

---

#### 5  Key Notebook Sections & Code Highlights

| Section                   | What it does                                                                                                   | Main functions                       |
| ------------------------- | -------------------------------------------------------------------------------------------------------------- | ------------------------------------ |
| **1. Imports & paths**    | Sets `DATA_DIR`, `RAW_DIR`, `EXTRACT_DIR` globals                                                              | `Pathlib`, `os`                      |
| **2. Untarring**          | Walks through `RAW_DIR`, opens each tarball and extracts to `EXTRACT_DIR/<GSE>/`                               | `tarfile.open(...).extractall()`     |
| **3. Decompression**      | Iterates over every `.gz` file inside each study directory; pipes to `shutil.copyfileobj` for streaming gunzip | `gzip.open`, `shutil`                |
| **4. Metadata ingestion** | Reads each `SraRunTable.csv`; re-indexes by `Sample Name` to ease joins                                        | `pd.read_csv`, `.set_index`          |
| **5. Matrix assembly**    | Reads text matrices into `pd.DataFrame`, aligns genes, stacks into a dict keyed by `GSE`                       | `pd.read_table`, `pd.concat(axis=1)` |
| **6. AnnData creation**   | Converts each combined DataFrame to `sc.AnnData`; attaches per-cell metadata                                   | `sc.AnnData(X)`, `.obs`, `.var`      |
| **7. Merge & QC**         | `anndata.concat` with `join='outer'`, computes simple QC fields (`n_genes`, `percent_mito`, `log_counts`)      | `adata.obs`, `.var`, `scanpy.pp.*`   |
| **8. Output**             | Persists final object (`combined_by_obs.h5ad`) and pickle of raw counts                                        | `adata.write_h5ad()`                 |

---

#### 6  Assumptions & Normalisation Choices

* **Gene identifiers** – all matrices store *Ensembl gene IDs*; symbols are looked up via `mygene` later in the pipeline.
* **Zero-filled outer join** – genes absent from a study are filled with zeros to preserve sparsity.
* **No doublet removal here** – initial QC only removes cells with <200 genes or >25 % mitochondrial reads; doublet detection deferred to the analysis notebook.
* **Raw counts retained** – no library size normalisation is applied during extraction; this ensures downstream Bayesian NB models receive integer counts with offset terms.

---

#### 7  Running the Notebook

1. **Activate** the project conda/env.

2. **Launch Jupyter** in `src/`:

   ```bash
   jupyter lab extract_files.ipynb
   ```

3. **Run all cells** (≈ 8–10 min on SSD with \~2 GB RAM peak).
   Progress prints every major step; OS-level errors (e.g., missing tar files) surface immediately.

4. **Verify outputs**:

   ```python
   import anndata as ad
   adata = ad.read_h5ad("../data/combined_by_obs.h5ad")
   adata
   # AnnData object with n_obs ≈ 96 000  |  n_vars ≈ 23 400
   ```

---

#### 8  Troubleshooting

| Symptom                     | Likely Cause                       | Fix                                                                 |
| --------------------------- | ---------------------------------- | ------------------------------------------------------------------- |
| `tarfile.ReadError`         | Corrupted or truncated GEO archive | Re-download the `*_RAW.tar` file                                    |
| Out-of-memory on extraction | Low RAM / huge study               | Edit `CHUNK_SIZE` constant to process counts in batches             |
| Gene mismatches after merge | Different annotation versions      | Supply a custom `ensembl_to_symbol.tsv` and run the re-mapping cell |

---

#### 9  Next Steps

* `src/bayesian_vi.ipynb` reads the merged AnnData, performs CPM normalisation, negative-binomial modelling in **PyMC**, and computes posterior expression trends.
* The saved gene length table is used as an offset in `pm.NegativeBinomial` to adjust for transcript length bias.

---

#### 10  Citation & Acknowledgements

* Original GEO datasets: see individual study DOIs in the project’s `docs/datasets.md`.
* Extraction logic adapted from the **Scanpy** tutorial (Wolf *et al.* Nature Methods 2018) and best practices outlined in **Orchestrating Single-Cell Analysis with Bioconductor** (Amezquita *et al.* 2020).
