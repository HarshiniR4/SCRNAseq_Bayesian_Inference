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

| Path / file               | Purpose                                         |
| ------------------------- | ----------------------------------------------- |
| **docs/**                 | Narrative docs & figures (see map below)        |
| **src/** (`bayes_scran/`) | Re-usable Python code, CLI entry-points         |
| **src/notebooks/**        | Exploratory notebooks (small subsets)           |
| **Snakefile** + `config/` | Reproducible Snakemake workflow                 |
| **environment.yml**       | Conda spec (PyMC 5, Scanpy 1.10, Snakemake ≥ 8) |
| **data/raw/** *(ignored)* | Place raw GSE98969 archives here                |
| **results/**              | Posterior draws, metrics, HTML dashboard        |
| **LICENSE**               | MIT                                             |

---

## Documentation map

| Doc file                            | Synopsis                                      |
| ----------------------------------- | --------------------------------------------- |
| **`docs/00_overview.md`**           | One-page problem statement & flowchart        |
| **`docs/01_data-prep.md`**          | Extraction, QC thresholds, Zarr chunking      |
| **`docs/02_model-spec.md`**         | Full maths (likelihood, priors, DAG)          |
| **`docs/03_inference-workflow.md`** | ADVI settings, GPU/CPU tips                   |
| **`docs/04_validation.md`**         | ELBO checks, posterior predictive diagnostics |
| **`docs/05_results.md`**            | Runtime, memory, key plots                    |
| **`docs/06_troubleshooting.md`**    | lib64 clashes, KaTeX, Snakemake FAQs          |

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
