# SCRNAseq\_Bayesian\_Inference

> **End‑to‑end Bayesian workflow for single‑cell RNA‑seq.**
> Harmonises six public microglia Alzheimer’s GEO studies and fits a hierarchical Negative‑Binomial model with minibatch variational inference.

---

## Table of Contents

1. [Project Overview](#project-overview)
2. [Data Flow](#data-flow)
3. [Repository Structure](#repository-structure)
4. [Installation](#installation)
5. [Quick Start](#quick-start)
6. [Detailed Pipeline](#detailed-pipeline)
7. [Data Management & Integrity](#data-management--integrity)
8. [Model Details](#model-details)
9. [Troubleshooting](#troubleshooting)
10. [Roadmap](#roadmap)
11. [Citation](#citation)
12. [License](#license)

---

## Project Overview

This repository assembles **six Alzheimer’s‑related GEO scRNA‑seq datasets** into a single `AnnData` object (≈186 k cells × 39 k genes) and applies a **hierarchical Negative‑Binomial Bayesian model** to discover condition‑specific expression programs.

Key design choices:

* **Streaming extraction** to keep <2 GB RAM during untar/gunzip.
* **Union‑gene merge** with sparse zero‑padding.
* **Minibatch ADVI** in **PyMC 5** for tractable Bayesian inference.
* **Conda‑locked** environments + checksum verification for full reproducibility.

---

## Data Flow

```mermaid
flowchart LR
    A[RAW GEO tarballs] --> B{extract_files.ipynb}
    B --> C[Per‑study .h5ad]
    C --> D{combine_adata.ipynb}
    D --> E[Merged AnnData (h5ad)]
    E --> F{neg_binom.ipynb}
    F --> G[Posterior trace + PPC]
```

* **Raw acquisition** → download `*_RAW.tar` + `SraRunTable.csv`.
* **Extraction** → stream untar │ gunzip → sparse matrices.
* **Merge** → outer gene union, QC, write `combined_by_obs.h5ad`.
* **Model fit** → hierarchical NB with library‑size & gene‑length offsets.

---

## Repository Structure

```text
SCRNAseq_Bayesian_Inference/
├── data/                # Git‑ignored artefacts
│   ├── raw/             # GEO tarballs (+ *.md5)
│   ├── extracted/       # Plain‑text matrices
│   ├── qc/              # PNGs + TSVs
│   └── combined_by_obs.h5ad
├── docs/                # Markdown docs (overview + notebook guides)
├── envs/                # Conda YAMLs (cpu.yml, gpu.yml, pymc.yml)
├── scripts/             # Helper CLI scripts
├── src/                 # Jupyter notebooks + utils
└── .gitignore
```

##  Repository contents

| Path / file                                                               | Purpose                                                    |
| ------------------------------------------------------------------------- | ---------------------------------------------------------- |
| **docs/Overview Documentation.md**                                        | One-page overview & high-level flowchart                   |
| **docs/cRNA seq- Datasets .md**                                           | Short descriptions of every GEO / SRA table                |
| **docs/Detailed Breakdown of *Bayesian Inference of Gene Expression*.md** | Full mathematical model **& DAG**                          |
| **docs/scRNASeq\_analysis\_documentation.docx**                           | Legacy notes (binary)                                      |
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

## Installation

```bash
# clone
git clone https://github.com/HarshiniR4/SCRNAseq_Bayesian_Inference.git
cd SCRNAseq_Bayesian_Inference

# CPU env (≈ 3 h model fit)
conda env create -f envs/cpu.yml
conda activate scrna_nb

# GPU env (CUDA 12 +, JAX)
conda env create -f envs/gpu.yml
conda activate scrna_nb
export JAX_PLATFORM_NAME=gpu
```
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

## Troubleshooting

| Symptom             | Fix                                                |
| ------------------- | -------------------------------------------------- |
| `tarfile.ReadError` | Re‑download corrupt archive.                       |
| Memory spike >8 GB  | Reduce HVG count or use backed mode.               |
| ELBO plateau        | Lower learning rate (`adam_lr=5e‑4`).              |
| CUDA not found      | Use `envs/cpu.yml` or install compatible CUDA ≥12. |

---

## Roadmap

* Zero‑Inflated NB variant.
* Batch‑aware random effects.
* DVC‑based data tracking.
* Gene‑set enrichment (fgsea) integration.

---

## Citation

If you use this code, please cite the original GEO datasets (see `docs/datasets.md`) and:

* Wolf et al., *SCANPY* (2018)
* Salvatier et al., *PyMC* (2016)
* Bürkner & Vehtari, Hierarchical NB for scRNA‑seq (2023)

---

## License

Licensed under the **MIT License**—see [`LICENSE`](LICENSE).

