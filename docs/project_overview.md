### **Project Overview – Pipeline**

> End-to-end pipeline for curating six Alzheimer’s-related GEO studies, harmonising them into a single **AnnData** object, and fitting a **hierarchical Negative-Binomial Bayesian model** with **minibatch variational inference**.
> 
>
> **Tech-stack pillars**
> *Data plumbing* → `tarfile`, `gzip`, **Scanpy/AnnData**, **pandas**
> 
> *Storage* → sparse **HDF5** (`.h5ad`) + on-disk backed mode
> 
> *Modelling* → **PyMC 5** + **Aesara/JAX**
> 
> *Reproducibility* → conda envs + lock files, deterministic random seeds, Git-ignored artefacts

---

## 1.  Data Lineage & Management Techniques

| Stage                                            | What Happens                                            | Key Techniques / Rationale                                                                                                                                                                                                                                                           |                                                                                                                        |
| ------------------------------------------------ | ------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------- |
| **Raw acquisition**                              | Download 6 × `*_RAW.tar` archives + GEO metadata tables | • `wget –c` with SHA-256 checksum verification  <br>• Optional **Aspera** script for faster transfer                                                                                                                                                                                 |                                                                                                                        |
| **Cold storage**                                 | Archives copied to `data/raw/` (Git-ignored)            | • Immutable filenames, date-stamped subfolder <br>• **ZFS snapshot** (if on NAS) for disaster recovery                                                                                                                                                                               |                                                                                                                        |
| **Archive streaming**<br>(`extract_files.ipynb`) | Untar → gunzip → plain-text matrices                    | • \`tarfile.open(..., mode="r                                                                                                                                                                                                                                                        | \*")`streams blocks, avoids full extraction to /tmp  <br>•`shutil.copyfileobj`+`gzip.open\` pipeline keeps RAM <200 MB |
| **Matrix parsing**                               | Convert text to CSR sparse matrices                     | • `pandas.read_table(chunksize=100_000)` + `scipy.sparse.csr_matrix`  <br>• Column ordering preserved exactly as read to retain barcode integrity                                                                                                                                    |                                                                                                                        |
| **Metadata ingestion**                           | Harmonise study & sample annotations                    | • Column standardisation map (`RAW_COLUMN_MAP`)  <br>• `CategoryDtype` applied early to cut memory                                                                                                                                                                                   |                                                                                                                        |
| **Gene ID harmonisation**                        | Ensembl IDs only                                        | • Lookup table created once (`ensembl_map.tsv`), hashed with MD5  <br>• Duplicate aliases resolved by priority list (Ensembl > HGNC > symbol)                                                                                                                                        |                                                                                                                        |
| **Per-study AnnData**                            | `.h5ad` with backed `r` mode                            | • Sparse `X` (≈98 % zeros)  <br>• `.raw` slot left blank to save disk                                                                                                                                                                                                                |                                                                                                                        |
| **Global merge**<br>(`combine_adata.ipynb`)      | Outer union of 39 241 genes                             | • **Two-pass strategy**: pass 1 reads only `.var_names` in backed mode to build master union; pass 2 loads each file, pads zeros, then concatenates via `anndata.concat(..., join='outer', merge='unique')`  <br>• Intermediate objects deleted + `gc.collect()` to cap RAM at ≈4 GB |                                                                                                                        |
| **Quality control**                              | Basic filters only (cells <200 genes, >25 % mito)       | • `scanpy.pp.calculate_qc_metrics`  <br>• QC thresholds logged to `logs/qc_summary.tsv`                                                                                                                                                                                              |                                                                                                                        |
| **Final AnnData**                                | `combined_by_obs.h5ad` (186 224 × 39 241)               | • Written with `compression="gzip", compression_opts=4` for 7 × disk reduction  <br>• Index MD5 stored as file attribute for integrity checks                                                                                                                                        |                                                                                                                        |

---

## 2.  Modelling & Computational Strategies

| Component                | Implementation Detail                                         | Why                                                               |
| ------------------------ | ------------------------------------------------------------- | ----------------------------------------------------------------- |
| **Library-size offset**  | `np.log1p(total_counts)` added to predictor                   | Normalises for sequencing depth without transforming counts       |
| **Gene-length offset**   | `np.log(gene_length / 1 kb)`                                  | Removes length bias (long genes accumulate more UMIs)             |
| **Hierarchical priors**  | Non-centred parameterisation for gene- and cell-level effects | Prevents funnel pathologies, speeds ADVI                          |
| **Dispersion**           | Gene-specific φ<sub>g</sub> \~ Half-Cauchy(2)                 | Captures over-dispersion typical of scRNA-seq                     |
| **Minibatching**         | 512 cells × 3 000 HVGs                                        | Fits in GPU/CPU cache (<1 GB) while still giving stable gradients |
| **ADVI**                 | 15 k updates (warm-up) → 60 k updates (main)                  | Convergence monitored by `CheckParametersConvergence(0.01)`       |
| **GPU acceleration**     | Optional: set `JAX_PLATFORM_NAME=gpu`                         | 4–8 × speed-up on A6000                                           |
| **Posterior predictive** | 100 replicates per cell                                       | Validates mean-variance trend & zero proportion                   |

---

## 3.  Repository Map (Expanded)

```
.
├── data/
│   ├── raw/                 # immutable GEO tarballs + checksums
│   ├── extracted/           # txt matrices (auto-cleanable)
│   ├── qc/                  # per-study QC PNGs + TSVs
│   └── combined_by_obs.h5ad # 6 GB (gzip-compressed HDF5)
│
├── docs/
│   ├── Overview Documentation.md   ← *this file*
│   ├── extract_files.md
│   ├── combine_adata.md
│   ├── neg_binom.md
│   └── datasets.md                # DOIs + accession metadata
│
├── envs/
│   ├── cpu.yml     # BLAS=openblas, pymc=5, no JAX
│   ├── gpu.yml     # cudatoolkit=12.4, jax[cuda], pymc=5
│   └── pymc.yml    # lean env used by CI
│
├── scripts/
│   ├── download_geo.sh            # wget/aspera + md5 check
│   └── run_pipeline.sh            # CLI wrapper for headless runs
│
├── src/
│   ├── extract_files.ipynb
│   ├── combine_adata.ipynb
│   ├── neg_binom.ipynb
│   ├── Bayesian Inference for Anndata.ipynb
│   └── utils/
│       ├── io_helpers.py          # streaming untar / gunzip
│       ├── gene_id_utils.py       # ID mapping + caching
│       └── minibatch_iter.py      # PyMC batch generator
└── .gitignore                     # excludes *.h5ad, *.tar, /data
```

---

## 4.  Execution Cheat-Sheet

```bash
# 0. Environment
conda env create -f envs/gpu.yml
conda activate scrna_nb
export JAX_PLATFORM_NAME=gpu               # set to 'cpu' if no GPU

# 1. Download raw data  (~20 GB)
bash scripts/download_geo.sh

# 2. Extract & preprocess  (~10 min, <2 GB RAM)
jupyter nbconvert --to notebook --execute src/extract_files.ipynb

# 3. Merge AnnData  (~5 min)
jupyter nbconvert --to notebook --execute src/combine_adata.ipynb

# 4. Fit model  (A6000: 30–40 min | CPU-only: 3–4 h)
jupyter nbconvert --to notebook --execute src/neg_binom.ipynb
```

All intermediate notebooks are executed **headless** and saved as `*_executed.ipynb` for provenance.

---

## 5.  Data-integrity & Reproducibility

* **Checksums** – Every binary artefact (\*.tar, \*.h5ad) carries an MD5 in a sidecar `.md5` file; `scripts/verify.sh` validates the full chain.
* **Deterministic seeds** – NumPy, Aesara, and JAX seeds are set at notebook start; the value is logged to `logs/runinfo.yaml`.
* **Environment lock** – `conda list --explicit > envs/lock.txt` after major changes.
* **CI smoke-test** – GitHub Action executes `run_pipeline.sh` on a small synthetic dataset each push.

---

## 6.  Planned Enhancements

| Milestone                    | Detail                                                  |
| ---------------------------- | ------------------------------------------------------- |
| **Zero-Inflated NB variant** | Compare WAIC / PSIS-LOO vs. current NB                  |
| **Batch-aware model**        | Add random intercept per sequencing run                 |
| **faiss-backed cache**       | Store gene-metadata lookups for instant re-runs         |
| **DVC pipeline**             | Track large artefacts + automatic re-execution triggers |

---
