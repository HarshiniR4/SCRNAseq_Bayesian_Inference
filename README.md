# SCRNAseq\_Bayesian\_Inference

> **End‑to‑end Bayesian workflow for single‑cell RNA‑seq.**
> Harmonises six public microglia Alzheimer’s GEO studies and fits a hierarchical Negative‑Binomial model with minibatch variational inference.


## Project Overview

This repository assembles **six Alzheimer’s‑related GEO scRNA‑seq datasets** into a single `AnnData` object (≈186 k cells × 39 k genes) and applies a **hierarchical Negative‑Binomial Bayesian model** to discover condition‑specific expression programs.

Key design choices:

* **Streaming extraction** to keep <2 GB RAM during untar/gunzip.
* **Union‑gene merge** with sparse zero‑padding.
* **Minibatch ADVI** in **PyMC 5** for tractable Bayesian inference.
* **Conda‑locked** environments + checksum verification for full reproducibility.

---

## Data Flow

```
RAW GEO tarballs → extract_files.ipynb → per‑study .h5ad
                                    ↓
                                    combine_adata.ipynb → combined_by_obs.h5ad
                                                            ↓
                                                            neg_binom.ipynb → posterior trace (.nc)
```

## Repository Structure

```
 SCRNAseq_Bayesian_Inference/
├─ data/                # large artefacts (git‑ignored)
│  ├─ raw/              #   GEO *.tar (immutable)
│  ├─ extracted/        #   temporary text matrices
│  ├─ combined_by_obs.h5ad
│  └─ nb_trace.nc        #   posterior draws (optional)
│
├─ docs/                # **human‑readable documentation**
│  ├─ Overview Documentation.md      # pipeline narrative
│  ├─ extract_files.md               # notebook‑level guide
│  ├─ combine_adata.md
│  ├─ neg_binom.md
│  └─ datasets.md                    # GEO accession table + DOIs
│
├─ envs/               # conda environments
│  ├─ cpu.yml
│  ├─ gpu.yml
│  └─ pymc.yml
│
├─ scripts/            # command‑line helpers (non‑notebook)
│  ├─ download_geo.sh           # wget/aspera + checksum
│  ├─ run_pipeline.sh           # headless execution chain
│  └─ verify_md5.sh             # integrity check for data/raw
│
├─ src/                # **executable notebooks & utils**
│  ├─ extract_files.ipynb
│  ├─ combine_adata.ipynb
│  ├─ neg_binom.ipynb
│  ├─ Bayesian Inference for Anndata.ipynb  # small SVI demo
│  └─ utils/
│      ├─ io_helpers.py          # streaming untar / gunzip
│      ├─ gene_id_utils.py       # Ensembl ↔ symbol mapping
│      └─ minibatch_iter.py      # PyMC minibatch generator
│
└─ README.md            # you are here
```

Key principles:

* **Separation of code vs. data** – `data/` is never version‑controlled; `.gitignore` enforces this.
* **One notebook → one markdown explainer** – every computational step is documented for auditability.
* **Utilities isolated** – pure‑Python helpers live in `src/utils/` to keep notebooks concise.

---

## Documentation Map

| Doc file                    | Scope                                             | Linked notebook       |
| --------------------------- | ------------------------------------------------- | --------------------- |
| `Overview Documentation.md` | Narrative of entire pipeline                      | –                     |
| `extract_files.md`          | Input expectations, step‑by‑step extraction logic | `extract_files.ipynb` |
| `combine_adata.md`          | Two‑pass merge strategy, zero‑padding details     | `combine_adata.ipynb` |
| `neg_binom.md`              | Full model spec, priors, VI settings              | `neg_binom.ipynb`     |
| `datasets.md`               | GEO accession metadata & citations                | none                  |

All docs live in **`docs/`** and are rendered natively by GitHub for quick browsing.

##  Repository contents

| Path / file                                                               | Purpose                                                    |
| ------------------------------------------------------------------------- | ---------------------------------------------------------- |
| **src/extract\_files.ipynb**                                              | Decompress `.tar`/`.gz`, build sparse CSR, attach metadata |
| **src/combine\_adata.ipynb**                                              | QC filters, gene-union realignment, write Zarr             |
| **src/Bayesian Inference for Anndata.ipynb**                              | ADVI on the entire dataset                                 |
| **src/VI\_Trace.ipynb**                                                   | ELBO trace & posterior diagnostics                         |
| **src/scRNAseq\_analysis.ipynb**                                          | Misc. exploration / plotting                               |
| **src/GSE* SraRunTable.csv*\*                                             | Sample-level metadata tables                               |
| **GSE123025\_Single\_myeloid\_1922\_cells\_processed\_data.csv.gz**       | Example processed matrix                                   |

*(Any large Zarr, `.h5ad`, or `.nc` artefacts should live in `results/` and be git-ignored.)*

---
## Mathematical model (mini version)

> See *Detailed Breakdown …Gene Expression\_.md* for full derivation and DAG.

### Likelihood

For cell *i* and gene *g*

$$
Y_{ig} \sim \mathrm{NegBinomial}\!\bigl(\mu_{ig},\phi_{g}\bigr),\qquad
\mu_{ig}=s_i\,e^{\eta_{ig}}.
$$

* $s_i$ — library-size factor
* $\phi_g$ — gene-specific inverse dispersion

### Linear predictor

$$
\begin{aligned}
\eta_{ig} &= \beta_{0g}
            + \beta_{1g}\,\mathrm{condition}_{i}
            + \gamma_{g}\,\log\!\bigl(\mathrm{geneLength}_{g}\bigr)
            + w_{ig}, \\[6pt]
w_{ig} &\sim \mathcal{N}\!\bigl(0,\sigma_w^{2}\bigr).
\end{aligned}
$$

### Priors

$$
\log s_i \sim \mathcal{N}\!\bigl(\mu_s,\sigma_s^{2}\bigr),\qquad
\beta_{\!\cdot g}\sim\mathcal{N}\!\bigl(0,\sigma_{\beta}^{2}I\bigr),\qquad
\phi_g\sim\mathrm{HalfCauchy}(\gamma).
$$

---

## DAG for Bayesian Inference 

![DAG](https://github.com/HarshiniR4/SCRNAseq_Bayesian_Inference/blob/main/images/Screenshot%202025-05-06%20120226.png)

--- 

## Detailed Pipeline

| Step               | Notebook                               | Highlights                                                                     |
| ------------------ | -------------------------------------- | ------------------------------------------------------------------------------ |
| **1 Extraction**   | `extract_files.ipynb`                  | Streaming untar + gzip; gene‑ID harmonisation; writes per‑study `.h5ad`.       |
| **2 Combine**      | `combine_adata.ipynb`                  | Two‑pass gene‑union; QC filters; 7× disk compression.                          |
| **3 Bayesian Fit** | `neg_binom.ipynb`                      | Minibatch ADVI (512 cells × 3 000 HVGs); hierarchical priors; GPU‑accelerated. |
| **4 SVI Demo**     | `Bayesian Inference for Anndata.ipynb` | Toy 1 k×500 subset; proof‑of‑concept SVI.                                      |

Full notebook‑level documentation lives in **`docs/`**.

---

## Data Management & Integrity

* **Checksums** – MD5 sidecar files for every binary; `scripts/verify.sh` validates chain.
* **Sparse HDF5** – `.h5ad` compressed with `gzip –4` (6 GB on disk).
* **Backed mode** – Large matrices read with `backed='r'` to cap RAM.
* **Deterministic seeds** – Logged to `logs/runinfo.yaml`.
* **CI smoke‑test** – GitHub Action runs pipeline on mock dataset each push.

---

## Model Details

* **Likelihood** – Negative‑Binomial with gene‑specific dispersion.
* **Offsets** – log(library size) + log(gene length / 1 kb).
* **Hierarchies** – non‑centred for gene (β₀, β₁) and cell (αᵢ) effects.
* **Inference** – PyMC 5 ADVI, `CheckParametersConvergence(0.01)`.
* **Diagnostics** – ELBO trace, posterior predictive check, ArviZ summary.
---

## License

Licensed under the **MIT License**—see [`LICENSE`](LICENSE).

