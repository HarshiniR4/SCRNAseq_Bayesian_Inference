# Bayesian Inference for Single-Cell RNA-Seq  
*A scalable PyMC 5 pipeline with Zarr storage and stochastic variational inference*

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

---

## Project snapshot

* **Data** 10× Chromium v3 counts (accession **GSE98969**): **186 224 cells × 39 241 genes** after QC  
* **Goal** Estimate gene-wise means & dispersions with full Bayesian uncertainty  
* **Model** Hierarchical Negative-Binomial with covariates and latent factors  
* **Inference engine** PyMC 5 + ADVI (mini-batch ≈ 128 cells × 128 genes)  
* **Peak RAM** ≈ 38 GB (CPU)  **Peak VRAM** 5 GB (A40)  
* **Run-time (full data)** ≈ 3 h 45 m on a single A40; ≈ 6 h 20 m CPU-only  
* **Outputs** Posterior draws (`*.nc`), QC report, HTML results dashboard

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

## Repository layout

SCRNAseq_Bayesian_Inference/
├── README.md
├── LICENSE
├── environment.yml
├── Snakefile
│
├── docs/
│   ├── 00_overview.md
│   ├── 01_datasets.md
│   ├── 02_data_prep.md
│   ├── 03_model_spec.md          ← DAG lives here
│   ├── 04_inference_workflow.md
│   ├── 05_validation.md
│   ├── 06_results.md
│   ├── 99_troubleshooting.md
│   └── figures/
│        ├── dag.svg
│        └── …
│
├── notebooks/
│   ├── preprocessing/
│   │   ├── extract_files.ipynb
│   │   └── combine_adata.ipynb
│   ├── inference/
│   │   ├── Bayesian_Inference_for_AnnData.ipynb
│   │   └── VI_trace_demo.ipynb
│   └── exploration/
│        └── scRNASeq_analysis.ipynb
│
├── src/
│   └── bayes_scran/
│        ├── __init__.py
│        ├── io_zarr.py
│        ├── qc.py
│        ├── model.py
│        └── cli.py
│
├── data/          # <ignored> raw MTX / FASTQ
└── results/       # <ignored> posterior draws, HTML reports

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

## Mathematical model

> Display math renders via GitHub KaTeX (public repos) — blank lines and **safe macros** only.

### Likelihood

For cell *i* and gene *g*

$$
Y_{ig} \sim \mathrm{NegBinomial}\!\bigl(\mu_{ig}, \phi_{g}\bigr),\qquad
\mu_{ig}=s_i\,\exp(\eta_{ig}).
$$

* $s_i$ — library-size factor
* $\phi_{g}$ — gene-specific inverse dispersion

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

### Hierarchical priors

$$
\log s_i \sim \mathcal{N}\!\bigl(\mu_s,\sigma_s^{2}\bigr),\qquad
\beta_{\!\cdot g}\sim\mathcal{N}\!\bigl(0,\sigma_{\beta}^{2} I\bigr),\qquad
\phi_{g} \sim \mathrm{HalfCauchy}(\gamma).
$$

*A full DAG appears in* `docs/Overview Documentation.md`.

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

## Contributing & license

* Pull requests welcome—new QC modules, inference engines, or docs fixes.
* Open issues for questions, bugs, or feature requests.

Licensed under the **MIT License**—see [`LICENSE`](LICENSE).

---
