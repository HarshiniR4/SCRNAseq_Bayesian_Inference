# Bayesian Inference for Single-Cell RNA-Seq  
*A scalable PyMC 5 pipeline with Zarr storage and stochastic variational inference*

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

---

## 1 Project at a glance

| Item | Value |
|------|-------|
| **Raw data** | GEO accession **GSE98969** (10× Chromium v3) |
| **Post-QC matrix** | 186 224 cells × 39 241 genes |
| **Model** | Hierarchical Negative-Binomial with covariates & latent factors |
| **Inference** | *Stochastic Variational Inference* (ADVI, PyMC 5) |
| **Peak RAM / VRAM** | 38 GB RAM · 5 GB (A40) |
| **Full-run wall-time** | ≈ 3 h 45 m (single A40) |
---

## Table of contents
1. [Quick start](#quick-start)  
2. [Repository layout](#repository-layout)  
3. [Documentation map](#documentation-map)  
4. [Mathematical model](#mathematical-model)  
5. [Results overview](#results-overview)  
6. [Contributing & license](#contributing--license)

---

## Quick start

```bash
# clone & set up
git clone https://github.com/HarshiniR4/SCRNAseq_Bayesian_Inference.git
cd SCRNAseq_Bayesian_Inference

conda env create -f environment.yml        # or: mamba env create …
conda activate scrnaseq-bayes
```
---

## 3 Repository contents

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

### Documentation map

| Doc file | Synopsis |
| -------- | -------- |
| **`docs/00_overview.md`** | Problem statement & high-level flowchart |
| **`docs/01_datasets.md`** | Description of raw GEO datasets & metadata |
| **`docs/02_data_prep.md`** | Extraction → QC → Zarr chunking |
| **`docs/03_model_spec.md`** | Full likelihood, priors & **DAG** |
| **`docs/04_inference_workflow.md`** | ADVI settings, PyTensor flags |
| **`docs/05_validation.md`** | ELBO, PPC, simulation recovery |
| **`docs/06_results.md`** | Runtime, memory, key posterior plots |
| **`docs/99_troubleshooting.md`** | lib64 clashes, KaTeX, Snakemake tips |

---


## 4 Mathematical model (mini version)

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

## Results overview

| Metric            | Value             |
| ----------------- | ----------------- |
| Cells × genes     | 186 224 × 39 241  |
| Peak RAM          | 38 GB             |
| Peak VRAM         | 5 GB (A40)        |
| ELBO plateau      | ≈ 18 k iterations |
| Held-out PPC RMSE | 0.14 ± 0.02       |

See `results/posterior.nc` (ArviZ-compatible) and the interactive **HTML dashboard** in `results/report.html`.

---

Licensed under the **MIT License**—see [`LICENSE`](LICENSE).

---
